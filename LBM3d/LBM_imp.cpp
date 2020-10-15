
#include "LBM3d.h"
#include <math.h>
#include <tbb/parallel_for.h>
//_______________________________________________________________________________________________________________________
///
/// \brief 3D vector
///
typedef struct
{
    real x;		//!< x component
    real y;		//!< y component
    real z;		//!< z component
}
vec3_t;

void   D3Q19( real* ex, real* ey, real* ez, int* oppos, real* wt )
{

    /// \brief Velocity vectors for D3Q19 as floating-point values
    const vec3_t vel3Dv[19] = {
      {  0,  0,  0 },		// zero direction
      { -1,  0,  0 },		// 6 directions with velocity 1
      {  1,  0,  0 },
      {  0, -1,  0 },
      {  0,  1,  0 },
      {  0,  0, -1 },
      {  0,  0,  1 },
      { -1, -1,  0 },		// 12 directions with velocity sqrt(2)
      { -1,  1,  0 },
      {  1, -1,  0 },
      {  1,  1,  0 },
      {  0, -1, -1 },
      {  0, -1,  1 },
      {  0,  1, -1 },
      {  0,  1,  1 },
      { -1,  0, -1 },
      { -1,  0,  1 },
      {  1,  0, -1 },
      {  1,  0,  1 },
    };

    const real init_weights3D[19] = { (real)1. / 3,
  (real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18,
  (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36 };



    for (int i = 0; i < 19; ++i) {
        ex[i] = vel3Dv[i].x;
        ey[i] = vel3Dv[i].y;
        ez[i] = vel3Dv[i].z;
        wt[i] = init_weights3D[i];
    }

 

    // define opposite (anti) aections (useful for implementing bounce back)

    const int invVel3D[19] = { 0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15 };
    for (int i = 0; i < 19; ++i) {
        oppos[i] = invVel3D[i];
    }
 
}









// this function updates the values of the distribution functions at all points along all directions
// carries out one lattice time-step (streaming + collision) in the algorithm

void collideAndStream(// READ-ONLY parameters (used by this function but not changed)
    const int NX, const int NY, const int NZ, const int Q, const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
    const real* ex, const real* ey, const real* ez, const int* oppos, const real* wt,
    // READ + WRITE parameters (get updated in this function)
    real* rho,         // density
    real* ux,         // X-velocity
    real* uy,         // Y-velocity
    real* uz,         // Z-velocity
    real* sigma,      // rate-of-strain
    real* f,          // distribution function
    real* feq,        // equilibrium distribution function
    real* f_new )      // new distribution function
{
    // loop over all interior voxels

    tbb::parallel_for( tbb::blocked_range<int>( 1, NX - 1 ),
        [&]( tbb::blocked_range<int> i_range ) {   for (int i = i_range.begin(); i < i_range.end(); ++i)


        //for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < NY - 1; j++)
        {
            for (int k = 1; k < NZ - 1; k++)
            {

                // natural index
                int index = i * NY * NZ + j * NZ + k;  // column-major ordering

                // calculate fluid viscosity based on the Reynolds number
                real kinematicViscosity = LID_VELOCITY * (real)NX / REYNOLDS_NUMBER;

                // calculate relaxation time tau
                real tau = 0.5 + 3.0 * kinematicViscosity;

                // collision
                for (int a = 0; a < Q; a++) {
                    int index_f = a + index * Q;
                    real edotu = ex[a] * ux[index] + ey[a] * uy[index] + ez[a] * uz[index];
                    real udotu = ux[index] * ux[index] + uy[index] * uy[index] + uz[index] * uz[index];
                    feq[index_f] = rho[index] * wt[a] * (1.0 + 3.0 * edotu + 4.5 * edotu * edotu - 1.5 * udotu);
                }

                // streaming from interior node points

                for (int a = 0; a < Q; a++) {

                    int index_f = a + index * Q;
                    int index_nbr = (i + ex[a]) * NY * NZ + (j + ey[a]) * NZ + (k + ez[a]);
                    int index_nbr_f = a + index_nbr * Q;
                    int indexoppos = oppos[a] + index * Q;

                    real tau_eff, tau_t, C_Smagorinsky;  // turbulence model parameters

                    C_Smagorinsky = 0.16;

                    // tau_t = additional contribution to the relaxation time 
                    //         because of the "eddy viscosity" model
                    // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                    // REFERENCE: Krafczyk M., Tolke J. and Luo L.-S. (2003)
                    //            Large-Eddy Simulations with a Multiple-Relaxation-Time LBE Model
                    //            International Journal of Modern Physics B, Vol.17, 33-39
                    // =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                    tau_t = 0.5 * (pow( pow( tau, 2.0 ) + 18.0 * pow( C_Smagorinsky, 2 ) * sigma[index], 0.5 ) - tau);

                    // the effective relaxation time accounts for the additional "eddy viscosity"
                    // effects. Note that tau_eff now varies from point to point in the domain, and is
                    // larger for large strain rates. If the strain rate is zero, tau_eff = 0 and we
                    // revert back to the original (laminar) LBM scheme where tau_eff = tau.

                    tau_eff = tau + tau_t;

                    // post-collision distribution at (i,j) along "a"
                    real f_plus = f[index_f] - (f[index_f] - feq[index_f]) / tau_eff;

                    int iS = i + ex[a];
                    int jS = j + ey[a];
                    int kS = k + ez[a];

                    if ((iS == 0) || (iS == NX - 1) || (jS == 0) || (jS == NY - 1) || (kS == 0) || (kS == NZ - 1)) {
                        // bounce back
                        real ubdote = ux[index_nbr] * ex[a] + uy[index_nbr] * ey[a] + uz[index_nbr] * ez[a];
                        f_new[indexoppos] = f_plus - 6.0 * DENSITY * wt[a] * ubdote;
                    }
                    else {
                        // stream to neighbor
                        f_new[index_nbr_f] = f_plus;
                    }
                }

            } // k
        } // j
    }
        }//i
    );
}

void macroVar( // READ-ONLY parameters (used by this function but not changed)
    const int NX, const int NY, const int NZ, const int Q, const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
    const real* ex, const real* ey, const real* ez, const int* oppos, const real* wt,
    // READ + WRITE parameters (get updated in this function)
    real* rho,         // density
    real* ux,         // X-velocity
    real* uy,         // Y-velocity
    real* uz,         // Z-velocity
    real* sigma,      // rate-of-strain
    real* f,          // distribution function
    real* feq,        // equilibrium distribution function
    real* f_new )      // new distribution function
{
    // loop over all interior voxels
    for (int i = 1; i < NX - 1; i++) {
        tbb::parallel_for( tbb::blocked_range<int>( 1, NY - 1 ),
            [&]( tbb::blocked_range<int> j_range ) {   for (int j = j_range.begin(); j < j_range.end(); ++j)
        {
            //for (int j = 1; j < N - 1; j++)

            for (int k = 1; k < NZ - 1; k++)
            {
                // natural index
                int index = i * NY * NZ + j * NZ + k;  // column-major ordering

                // push f_new into f
                for (int a = 0; a < Q; a++) {
                    int index_f = a + index * Q;
                    f[index_f] = f_new[index_f];
                }

                // update density at interior nodes
                rho[index] = 0.0;
                for (int a = 0; a < Q; a++) {
                    int index_f = a + index * Q;
                    rho[index] += f_new[index_f];
                }

                // update velocity at interior nodes
                real velx = 0.0;
                real vely = 0.0;
                real velz = 0.0;
                for (int a = 0; a < Q; a++) {
                    int index_f = a + index * Q;
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
                for (int a = 1; a < Q; a++)
                {
                    int index_f = a + index * Q;

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
                    + pow( sum_yx, 2 ) + pow( sum_yy, 2 ) + pow( sum_yz, 2 )
                    + pow( sum_zx, 2 ) + pow( sum_zy, 2 ) + pow( sum_zz, 2 );

                sigma[index] = pow( sigma[index], 0.5 );

            } //k
        }//j
            } );


    }
}




void initialize( const int NX, const int NY, const int NZ, const int Q, const real DENSITY, const real LID_VELOCITY,
    real* ex, real* ey, real* ez, int* oppos, real* wt,
    real* rho, real* ux, real* uy, real* uz, real* sigma,
    real* f, real* feq, real* f_new )
{
    // loop over all voxels
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {

                // natural index for location (i,j)

                int index = i * NY * NZ + j * NZ + k;  // column-ordering

                // initialize density and velocity fields inside the cavity

                rho[index] = DENSITY;   // density
                ux[index] = 0.0;       // x-component of velocity
                uy[index] = 0.0;       // y-component of velocity
                uz[index] = 0.0;       // z-component of velocity
                sigma[index] = 0.0;       // rate-of-strain field

                // specify boundary condition for the moving lid

                //if (j == N - 1) ux[index] = LID_VELOCITY;
               
                if ((j == NY-1)   ) { ux[index] = LID_VELOCITY; }

                // assign initial values for distribution functions
                // along various aections using equilibriu, functions

                for (int a = 0; a < Q; a++) {

                    int index_f = a + index * Q;

                    real edotu = ex[a] * ux[index] + ey[a] * uy[index] + ez[a] * uz[index];
                    real udotu = ux[index] * ux[index] + uy[index] * uy[index] + uz[index] * uz[index];
                     
                    feq[index_f] = rho[index] * wt[a] * (1.0 + 3.0 * edotu + 4.5 * edotu * edotu - 1.5 * udotu);
                    f[index_f] = feq[index_f];
                    f_new[index_f] = feq[index_f];

                }

            }
        }
    }
}