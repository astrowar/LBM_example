#include "cuda_runtime.h"

#include "lbm_gpu.h"
#include <math.h>
 
#include <stdio.h>

#define Qcc (19)

#define THREADS_PER_BLOCK 512


#define TILEDIM_X 2
#define TILEDIM_Y 8
#define TILEDIM_Z 8

#define CELLDIM_X (TILEDIM_X+2)
#define CELLDIM_Y (TILEDIM_Y+2)
#define CELLDIM_Z (TILEDIM_Z+2)


__device__ void collideAndStream_cuda_ijk( const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
	int i, int j, int k, size_t nx, size_t ny, size_t nz,
	 
	real* rho,         // density
	const real* ux,         // X-velocity
	const real* uy,         // Y-velocity
	const real* uz,         // Z-velocity
	real* sigma,      // rate-of-strain
	const real* f,          // distribution function
	real* feq,        // equilibrium distribution function
	real* f_new )      // new distribution function

{

	real ex[19] = { 0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000  };
	real ey[19] = { 0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000  };
	real ez[19] = { 0.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000  };
	real wt[19] = { 0.333333 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778  };
	int oppos[19] = { 0 ,2 ,1 ,4 ,3 ,6 ,5 ,10 ,9 ,8 ,7 ,14 ,13 ,12 ,11 ,18 ,17 ,16 ,15   };


	// natural index
	int index = i * ny * nz + j * nz + k;  // column-major ordering

	// calculate fluid viscosity based on the Reynolds number
	real kinematicViscosity = LID_VELOCITY * (real)nx / REYNOLDS_NUMBER;

	// calculate relaxation time tau
	real tau = 0.5 + 3.0 * kinematicViscosity;
	 real udotu = ux[index] * ux[index] + uy[index] * uy[index] + uz[index] * uz[index];
	 real roi = rho[index];
	// collision
	for (int a = 0; a < Qcc; a++) {
		int index_f = a + index * Qcc;
		real edotu = ex[a] * ux[index] + ey[a] * uy[index] + ez[a] * uz[index];		
		feq[index_f] = roi  * wt[a] * (1.0f + 3.0f * edotu + 4.5 * edotu * edotu - 1.5 * udotu);
	}

	// streaming from interior node points
 

	// tau_t = additional contribution to the relaxation time 
	//         because of the "eddy viscosity" model
	// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	// REFERENCE: Krafczyk M., Tolke J. and Luo L.-S. (2003)
	//            Large-Eddy Simulations with a Multiple-Relaxation-Time LBE Model
	//            International Journal of Modern Physics B, Vol.17, 33-39
	// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	real C_Smagorinsky = 0.16;
	 real tau_t = 0.5 * (powf( powf( tau, 2.0 ) + 18.0 * powf( C_Smagorinsky, 2 ) * sigma[index], 0.5 ) - tau);
	for (int a = 0; a < Qcc; a++) {

		int index_f = a + index * Qcc;
		int index_nbr = (i + ex[a]) * ny * nz + (j + ey[a]) * nz + (k + ez[a]);
		int index_nbr_f = a + index_nbr * Qcc;
		int indexoppos = oppos[a] + index * Qcc;

		real tau_eff; // turbulence model parameters

		

		// tau_t = additional contribution to the relaxation time 
		//         because of the "eddy viscosity" model
		// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
		// REFERENCE: Krafczyk M., Tolke J. and Luo L.-S. (2003)
		//            Large-Eddy Simulations with a Multiple-Relaxation-Time LBE Model
		//            International Journal of Modern Physics B, Vol.17, 33-39
		// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
		
		

		// the effective relaxation time accounts for the additional "eddy viscosity"
		// effects. Note that tau_eff now varies from point to point in the domain, and is
		// larger for large strain rates. If the strain rate is zero, tau_eff = 0 and we
		// revert back to the original (laminar) LBM scheme where tau_eff = tau.

		tau_eff = tau + tau_t;

		// post-collision distribution at (i,j) along "a"
		real f_plus = f[index_f] -  (f[index_f] - feq[index_f]) / tau_eff;

		//f_plus = f_plus * 2.0;

		int iS = i + ex[a];
		int jS = j + ey[a];
		int kS = k + ez[a];

		if ((iS == 0) || (iS == nx - 1) || (jS == 0) || (jS == ny - 1) || (kS == 0) || (kS == nz - 1)) {
			// bounce back
			real ubdote = ux[index_nbr] * ex[a] + uy[index_nbr] * ey[a] + uz[index_nbr] * ez[a];
			f_new[indexoppos] = f_plus - 6.0 * DENSITY * wt[a] * ubdote;
		}
		else {
			// stream to neighbor
			f_new[index_nbr_f] = f_plus;
		}
	}


	 
}





__device__ void collideAndStream_cuda_ijk_local( const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
	int i, int j, int k, size_t nx, size_t ny, size_t nz,
	const real* ex, const real* ey, const real* ez, const int* oppos, const real* wt,
	real* rho,         // density
	real* ux,         // X-velocity
	real* uy,         // Y-velocity
	real* uz,         // Z-velocity
	real* sigma,      // rate-of-strain
	real* f,          // distribution function
	real* feq,        // equilibrium distribution function
	real* f_new_local )      // new distribution function

{
	int i_local = threadIdx.x + 1;
	int j_local = threadIdx.y + 1;
	int k_local = threadIdx.z + 1;
	int  a = 0;
	// natural index
	int index = i * ny * nz + j * nz + k;  // column-major ordering
	int index_local = i_local * (CELLDIM_Y * CELLDIM_Z) + j_local * CELLDIM_Z + k_local;  // column-major ordering


	// calculate fluid viscosity based on the Reynolds number
	real kinematicViscosity = LID_VELOCITY * (real)nx / REYNOLDS_NUMBER;

	// calculate relaxation time tau
	real tau = 0.5 + 3.0 * kinematicViscosity;

	// collision
	for (  a = 0; a < Qcc; a++) {
		int index_f = a + index * Qcc;
		real edotu = ex[a] * ux[index] + ey[a] * uy[index] + ez[a] * uz[index];
		real udotu = ux[index] * ux[index] + uy[index] * uy[index] + uz[index] * uz[index];
		feq[index_f] = rho[index] * wt[a] * (1.0f + 3.0f * edotu + 4.5 * edotu * edotu - 1.5 * udotu);
	}

	// streaming from interior node points

	for (a = 0; a < Qcc; a++) {

		int index_f = a + index * Qcc;
		int index_nbr = (i + ex[a]) * ny * nz + (j + ey[a]) * nz + (k + ez[a]);
		int index_nbr_f = a + index_nbr * Qcc;
		int indexoppos = oppos[a] + index * Qcc;


		int index_local_nbr = (i_local + int( ex[a] )) * (CELLDIM_Y * CELLDIM_Z) + (j_local + int( ey[a] )) * (CELLDIM_Z)+(k_local + int( ez[a] ));
		int index_local_nbr_f = a + index_local_nbr * Qcc;
		int indexoppos_local = oppos[a] + index_local * Qcc;

		real tau_eff, tau_t, C_Smagorinsky;  // turbulence model parameters
		C_Smagorinsky = 0.16;
		tau_t = 0.5 * (powf( powf( tau, 2.0 ) + 18.0 * powf( C_Smagorinsky, 2 ) * sigma[index], 0.5 ) - tau);
		tau_eff = tau + tau_t;
		// post-collision distribution at (i,j) along "a"
		real f_plus = f[index_f] - (f[index_f] - feq[index_f]) / tau_eff;
		int iS = i + int( ex[a] );
		int jS = j + int( ey[a] );
		int kS = k + int( ez[a] );
		f_new_local[index_local_nbr_f] = 0;
		if ((iS == 0) || (iS == nx - 1) || (jS == 0) || (jS == ny - 1) || (kS == 0) || (kS == nz - 1)) {
			// bounce back
			real ubdote = ux[index_nbr] * ex[a] + uy[index_nbr] * ey[a] + uz[index_nbr] * ez[a];
			f_new_local[indexoppos_local] = f_plus - 6.0f * DENSITY * wt[a] * ubdote;
		}
		else {
			// stream to neighbor
			 f_new_local[index_local_nbr_f] = f_plus;
		}

	}
}




 

	__global__  void   collideAndStream_cuda(
		const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
		size_t nx, size_t ny, size_t nz,
		const real * ex, const real * ey, const real * ez, const int* oppos, const real * wt,
		real * rho,         // density
		real * ux,         // X-velocity
		real * uy,         // Y-velocity
		real * uz,         // Z-velocity
		real * sigma,      // rate-of-strain
		real * f,          // distribution function
		real * feq,        // equilibrium distribution function
		real * f_new
	) {
		int index = blockIdx.x * blockDim.x + threadIdx.x;
		int stridex = blockDim.x * gridDim.x;

		int indey = blockIdx.y * blockDim.y + threadIdx.y;
		int stridey = blockDim.y * gridDim.y;

		int indez = blockIdx.z * blockDim.z + threadIdx.z;
		int stridez = blockDim.z * gridDim.z;

		for (int i = index; i < nx - 1; i += stridex)
			for (int j = indey; j < ny - 1; j += stridey)
				for (int k = indez; k < nz - 1; k += stridez)
				{
					if ((i > 0) && (j > 0) && (k > 0))
					
						collideAndStream_cuda_ijk( DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, i, j, k, nx, ny, nz, rho, ux, uy, uz, sigma, f, feq, f_new );
				}
	}

	// this function updates the values of the distribution functions at all points along all directions
	// carries out one lattice time-step (streaming + collision) in the algorithm

	void  LBMGrid::collideAndStream(// READ-ONLY parameters (used by this function but not changed)
		const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER
	)
		// new distribution function
	{

		dim3 grid( nx / TILEDIM_X + 1, ny / TILEDIM_Y + 1, nz / TILEDIM_Z + 1 );
		dim3 block( TILEDIM_X, TILEDIM_Y, TILEDIM_Z );
		collideAndStream_cuda << < grid, block >> > (DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, nx, ny, nz, ex, ey, ez, oppos, wt, rho, ux, uy, uz, sigma, f, feq, f_new);
		cudaDeviceSynchronize();
	}

