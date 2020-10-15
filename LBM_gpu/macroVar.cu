#include "cuda_runtime.h"

#include "lbm_gpu.h"
#include <math.h>

#define Qcc (19)

#define THREADS_PER_BLOCK 512
#define TILEDIM_X 4
#define TILEDIM_Y 8
#define TILEDIM_Z 8

#define CELLDIM_X (TILEDIM_X+2)
#define CELLDIM_Y (TILEDIM_Y+2)
#define CELLDIM_Z (TILEDIM_Z+2)

__device__ void macroVar_cuda_ijk( const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
	int i, int j, int k, size_t nx, size_t ny, size_t nz, 
	 
	real* rho,         // density
	real* ux,         // X-velocity
	real* uy,         // Y-velocity
	real* uz,         // Z-velocity
	real* sigma,      // rate-of-strain
	real* f,          // distribution function
	real* feq,        // equilibrium distribution function
	real* f_new )      // new distribution function

{

	real ex[19] = { 0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 };
	real ey[19] = { 0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 };
	real ez[19] = { 0.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 };
	real wt[19] = { 0.333333 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 };
	int oppos[19] = { 0 ,2 ,1 ,4 ,3 ,6 ,5 ,10 ,9 ,8 ,7 ,14 ,13 ,12 ,11 ,18 ,17 ,16 ,15   };


	int index = i * ny * nz + j * nz + k;  // column-major ordering


			  // push f_new into f
	for (int a = 0; a < Qcc; a++) {
		int index_f = a + index * Qcc;
		f[index_f] = f_new[index_f];
	}

	// update density at interior nodes
	rho[index] = 0.0;
	for (int a = 0; a < Qcc; a++) {
		int index_f = a + index * Qcc;
		rho[index] += f_new[index_f];
	}

	// update velocity at interior nodes
	real velx = 0.0;
	real vely = 0.0;
	real velz = 0.0;
	for (int a = 0; a < Qcc; a++) {
		int index_f = a + index * Qcc;
		velx += f_new[index_f] * ex[a];
		vely += f_new[index_f] * ey[a];
		velz += f_new[index_f] * ez[a];
	}
	ux[index] = velx / rho[index];
	uy[index] = vely / rho[index];
	uz[index] = velz / rho[index];

	// update the rate-of-strain field
	real sum_xx = 0.0, sum_xy = 0.0, sum_xz = 0.0;
	real sum_yx = 0.0, sum_yy = 0.0, sum_yz = 0.0;
	real sum_zx = 0.0, sum_zy = 0.0, sum_zz = 0.0;
	for (int a = 1; a < Qcc; a++)
	{
		int index_f = a + index * Qcc;

		real dFF = (f_new[index_f] - feq[index_f]);
		sum_xx = sum_xx + dFF * ex[a] * ex[a];
		sum_xy = sum_xy + dFF * ex[a] * ey[a];
		sum_xz = sum_xz + dFF * ex[a] * ez[a];
		sum_yx = sum_yx + dFF * ey[a] * ex[a];;
		sum_yy = sum_yy + dFF * ey[a] * ey[a];
		sum_yz = sum_yz + dFF * ey[a] * ez[a];
		sum_zx = sum_zx + dFF * ez[a] * ex[a];
		sum_zy = sum_zy + dFF * ez[a] * ey[a];
		sum_zz = sum_zz + dFF * ez[a] * ez[a];
	}

	// evaluate |S| (magnitude of the strain-rate)
	sigma[index] = pow( sum_xx, 2 ) + pow( sum_xy, 2 ) + pow( sum_xz, 2 )
		+ powf( sum_yx, 2 ) + powf( sum_yy, 2 ) + powf( sum_yz, 2 )
		+ powf( sum_zx, 2 ) + powf( sum_zy, 2 ) + powf( sum_zz, 2 );

	sigma[index] = powf( sigma[index], 0.5 );

}


//LOCAL SHARED


__device__ void macroVar_cuda_ijk_local( const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
	int i, int j, int k, size_t nx, size_t ny, size_t nz,

 
	real* rho_local,         // density
	real* ux,         // X-velocity
	real* uy,         // Y-velocity
	real* uz,         // Z-velocity
	real* sigma,      // rate-of-strain
	real* f_local,          // distribution function
	real* feq_local,        // equilibrium distribution function
	real* f_new_local )      // new distribution function

{

	real ex[19] = { 0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 };
	real ey[19] = { 0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,-1.000000 ,1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 };
	real ez[19] = { 0.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 ,-1.000000 ,1.000000 };
	real wt[19] = { 0.333333 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.055556 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 ,0.027778 };
	int oppos[19] = { 0 ,2 ,1 ,4 ,3 ,6 ,5 ,10 ,9 ,8 ,7 ,14 ,13 ,12 ,11 ,18 ,17 ,16 ,15   };



	int i_local = threadIdx.x + 1;
	int j_local = threadIdx.y + 1;
	int k_local = threadIdx.z + 1;

	int index = i * ny * nz + j * nz + k;  // column-major ordering

	int index_local = i_local * (CELLDIM_Y * CELLDIM_Z) + j_local * CELLDIM_Z + k_local;  // column-major ordering
	//int index_local = (threadIdx.x + 1) * CELLDIM_Y * CELLDIM_Z + (threadIdx.y + 1) * CELLDIM_Z + (threadIdx.z + 1);  // column-major ordering

			  // push f_new into f
	for (int a = 0; a < Qcc; a++) {
		int index_local_f = a + index_local * Qcc;
		f_local[index_local_f] = f_new_local[index_local_f];
	}

	// update density at interior nodes
	rho_local[index_local] = 0.0;
	for (int a = 0; a < Qcc; a++) {
		int index_local_f = a + index_local * Qcc;
		rho_local[index_local] +=  f_new_local[index_local_f];
	}

	// update velocity at interior nodes
	real velx = 0.0;
	real vely = 0.0;
	real velz = 0.0;
	for (int a = 0; a < Qcc; a++) {
		int index_local_f = a + index_local * Qcc;
		velx += f_new_local[index_local_f] * ex[a];
		vely += f_new_local[index_local_f] * ey[a];
		velz += f_new_local[index_local_f] * ez[a];
	}
	ux[index] = velx / rho_local[index_local];
	uy[index] = vely / rho_local[index_local];
	uz[index] = velz / rho_local[index_local];

	// update the rate-of-strain field
	real sum_xx = 0.0, sum_xy = 0.0, sum_xz = 0.0;
	real sum_yx = 0.0, sum_yy = 0.0, sum_yz = 0.0;
	real sum_zx = 0.0, sum_zy = 0.0, sum_zz = 0.0;
	for (int a = 1; a < Qcc; a++)
	{
		int index_local_f = a + index_local * Qcc;

		real dFF = (f_new_local[index_local_f] - feq_local[index_local_f]);
		sum_xx = sum_xx + dFF * ex[a] * ex[a];
		sum_xy = sum_xy + dFF * ex[a] * ey[a];
		sum_xz = sum_xz + dFF * ex[a] * ez[a];
		sum_yx = sum_yx + dFF * ey[a] * ex[a];;
		sum_yy = sum_yy + dFF * ey[a] * ey[a];
		sum_yz = sum_yz + dFF * ey[a] * ez[a];
		sum_zx = sum_zx + dFF * ez[a] * ex[a];
		sum_zy = sum_zy + dFF * ez[a] * ey[a];
		sum_zz = sum_zz + dFF * ez[a] * ez[a];
	}

	// evaluate |S| (magnitude of the strain-rate)
	sigma[index] = powf( sum_xx, 2 ) + powf( sum_xy, 2 ) + powf( sum_xz, 2 )
		+ powf( sum_yx, 2 ) + powf( sum_yy, 2 ) + powf( sum_yz, 2 )
		+ powf( sum_zx, 2 ) + powf( sum_zy, 2 ) + powf( sum_zz, 2 );

	sigma[index] = powf( sigma[index], 0.5 );

}


//
//__global__  void   macroVar_cuda_local(
//	const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
//	size_t nx, size_t ny, size_t nz,
//	const real* ex, const real* ey, const real* ez, const int* oppos, const real* wt,
//	real* rho,         // density
//	real* ux,         // X-velocity
//	real* uy,         // Y-velocity
//	real* uz,         // Z-velocity
//	real* sigma,      // rate-of-strain
//	real* f,          // distribution function
//	real* feq,        // equilibrium distribution function
//	real* f_new
//) {
//
//	int inde_x = blockIdx.x * blockDim.x + threadIdx.x;
//	int stridex = blockDim.x * gridDim.x;
//
//	int inde_y = blockIdx.y * blockDim.y + threadIdx.y;
//	int stridey = blockDim.y * gridDim.y;
//
//	int inde_z = blockIdx.z * blockDim.z + threadIdx.z;
//	int stridez = blockDim.z * gridDim.z;
//	__shared__ float f_local[CELLDIM_X * CELLDIM_Y * CELLDIM_Z * Qcc];
//	__shared__ float feq_local[CELLDIM_X * CELLDIM_Y * CELLDIM_Z * Qcc];
//	__shared__ float f_new_local[CELLDIM_X * CELLDIM_Y * CELLDIM_Z * Qcc];
//	__shared__ float rho_local[CELLDIM_X * CELLDIM_Y * CELLDIM_Z];
//
//
//
//	//copy from global to local
//
//
//
//	for (int i = inde_x; i < nx - 1; i += stridex)
//		for (int j = inde_y; j < ny - 1; j += stridey)
//			for (int k = inde_z; k < nz - 1; k += stridez)
//			{
//				//copy from global to local
//				int index = i * ny * nz + j * nz + k;  // column-major ordering
//				int index_local = (threadIdx.x + 1) * CELLDIM_Y * CELLDIM_Z + (threadIdx.y + 1) * CELLDIM_Z + (threadIdx.z + 1);  // column-major ordering
//				for (int a = 0; a < Qcc; a++) {
//					int index_f = a + index * Qcc;
//					int index_f_local = a + index_local * Qcc;
//					f_local[index_f_local] = f[index_f];
//					feq_local[index_f_local] = feq[index_f];
//					f_new_local[index_f_local] = f_new[index_f];
//				}
//				//ux_local[index_local] = ux[index];
//			   // uy_local[index_local] = uy[index];
//			   // uz_local[index_local] = uz[index];
//				//sigma_local[index_local] = sigma[index];
//				rho_local[index_local] = rho[index];
//
//				//copy the boundaries
//				for (int di = -1; di <= 1; di += 1)
//					for (int dj = -1; dj <= 1; dj += 1)
//						for (int dk = -1; dk <= 1; dk += 1)
//						{
//							int thx = threadIdx.x + 1 + di;
//							int thy = threadIdx.y + 1 + dj;
//							int thz = threadIdx.z + 1 + dk;
//
//							if ((thx == 0) || (thy == 0) || (thz == 0) || (thx == CELLDIM_X - 1) || (thy == CELLDIM_Y - 1) || (thz == CELLDIM_Z - 1))
//								if ((i + di >= 0) && (j + dj >= 0) && (k + dk >= 0) && (i + di < nx) && (j + dj < ny) && (k + dk < nz))
//								{
//									int index = (i + di) * ny * nz + (j + dj) * nz + (k + dk);
//									int index_local = (thx)*CELLDIM_Y * CELLDIM_Z + (thy)*CELLDIM_Z + (thz);  // column-major ordering
//
//									for (int a = 0; a < Qcc; a++) {
//										int index_f = a + index * Qcc;
//										int index_f_local = a + index_local * Qcc;
//										f_local[index_f_local] = f[index_f];
//										feq_local[index_f_local] = feq[index_f];
//										f_new_local[index_f_local] = f_new[index_f];
//									}
//									rho_local[index_local] = rho[index];
//								}
//
//						}
//
//
//				__syncthreads();
//
//				if ((i > 0) && (j > 0) && (k > 0))
//					macroVar_cuda_ijk_local( DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, i, j, k, nx, ny, nz, 
//						rho_local, ux, uy, uz, sigma, f_local, feq_local, f_new_local );
//
//				__syncthreads();
//
//				//copy back to global memory
//				for (int a = 0; a < Qcc; a++) {
//					int index_f = a + index * Qcc;
//					int index_f_local = a + index_local * Qcc;
//					f[index_f] = f_local[index_f_local];
//					feq[index_f] = feq_local[index_f_local];
//					f_new[index_f] = f_new_local[index_f_local];
//				}
//				//ux[index] = ux_local[index_local] ;
//				//uy[index] =  uy_local[index_local]  ;
//				//uz[index] = uz_local[index_local]  ;
//				//sigma[index] = sigma_local[index_local] ;
//				rho[index] = rho_local[index_local];
//
//			}
//}


__global__  void   macroVar_cuda_raw(
	const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
	size_t nx, size_t ny, size_t nz,
	const real* ex, const real* ey, const real* ez, const int* oppos, const real* wt,
	real* rho,         // density
	real* ux,         // X-velocity
	real* uy,         // Y-velocity
	real* uz,         // Z-velocity
	real* sigma,      // rate-of-strain
	real* f,          // distribution function
	real* feq,        // equilibrium distribution function
	real* f_new
) {

	int inde_x = blockIdx.x * blockDim.x + threadIdx.x;
	int stridex = blockDim.x * gridDim.x;

	int inde_y = blockIdx.y * blockDim.y + threadIdx.y;
	int stridey = blockDim.y * gridDim.y;

	int inde_z = blockIdx.z * blockDim.z + threadIdx.z;
	int stridez = blockDim.z * gridDim.z;
 



	//copy from global to local



	for (int i = inde_x; i < nx - 1; i += stridex)
		for (int j = inde_y; j < ny - 1; j += stridey)
			for (int k = inde_z; k < nz - 1; k += stridez)
			{
  
				if ((i > 0) && (j > 0) && (k > 0))
					macroVar_cuda_ijk ( DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, i, j, k, nx, ny, nz,
						rho , ux, uy, uz, sigma, f , feq , f_new  ); 
			}
}


void  LBMGrid::macroVar(// READ-ONLY parameters (used by this function but not changed)
	const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER, bool copy
)
// new distribution function
{
	dim3 grid( nx / TILEDIM_X + 1, ny / TILEDIM_Y + 1, nz / TILEDIM_Z + 1 );
	dim3 block( TILEDIM_X, TILEDIM_Y, TILEDIM_Z );
	//macroVar_cuda << < grid, block >> > (DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, nx, ny, nz, ex, ey, ez, oppos, wt, rho, ux, uy, uz, sigma, f, feq, f_new);
	macroVar_cuda_raw << < grid, block >> > (DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, nx, ny, nz, ex, ey, ez, oppos, wt, rho, ux, uy, uz, sigma, f, feq, f_new);
	

	if (copy) {
		cudaDeviceSynchronize();
		cudaMemcpy( cpu_ux, ux, sizeof( real ) * nx * ny * nz, cudaMemcpyDeviceToHost );
		cudaMemcpy( cpu_uy, uy, sizeof( real ) * nx * ny * nz, cudaMemcpyDeviceToHost );
		cudaMemcpy( cpu_uz, uz, sizeof( real ) * nx * ny * nz, cudaMemcpyDeviceToHost );
	}
}
