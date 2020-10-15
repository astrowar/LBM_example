#include "cuda_runtime.h"
#include "lbm_gpu.h"

#define Qcc (19)
void LBMGrid::initialize( real DENSITY, real  LID_VELOCITY )
{
    // distribution functions
    real* _f = new real[nx* ny *nz * Qcc];
    real* _feq = new real[nx * ny * nz * Qcc];
    real* _f_new = new real[nx * ny * nz * Qcc];

    // density and velocity
    real* _rho = new real[nx * ny * nz];
    real* _ux = new real[nx * ny * nz];
    real* _uy = new real[nx * ny * nz];
    real* _uz = new real[nx * ny * nz];

    real* _ex = new real[Qcc];
    real* _ey = new real[Qcc];
    real* _ez = new real[Qcc];

    real* _wt = new real[Qcc];


    // rate-of-strain
    real* _sigma = new real[nx * ny * nz];

    _ex[0] = 0.000000;
    _ey[0] = 0.000000;
    _ez[0] = 0.000000;
    _ex[1] = -1.000000;
    _ey[1] = 0.000000;
    _ez[1] = 0.000000;
    _ex[2] = 1.000000;
    _ey[2] = 0.000000;
    _ez[2] = 0.000000;
    _ex[3] = 0.000000;
    _ey[3] = -1.000000;
    _ez[3] = 0.000000;
    _ex[4] = 0.000000;
    _ey[4] = 1.000000;
    _ez[4] = 0.000000;
    _ex[5] = 0.000000;
    _ey[5] = 0.000000;
    _ez[5] = -1.000000;
    _ex[6] = 0.000000;
    _ey[6] = 0.000000;
    _ez[6] = 1.000000;
    _ex[7] = -1.000000;
    _ey[7] = -1.000000;
    _ez[7] = 0.000000;
    _ex[8] = -1.000000;
    _ey[8] = 1.000000;
    _ez[8] = 0.000000;
    _ex[9] = 1.000000;
    _ey[9] = -1.000000;
    _ez[9] = 0.000000;
    _ex[10] = 1.000000;
    _ey[10] = 1.000000;
    _ez[10] = 0.000000;
    _ex[11] = 0.000000;
    _ey[11] = -1.000000;
    _ez[11] = -1.000000;
    _ex[12] = 0.000000;
    _ey[12] = -1.000000;
    _ez[12] = 1.000000;
    _ex[13] = 0.000000;
    _ey[13] = 1.000000;
    _ez[13] = -1.000000;
    _ex[14] = 0.000000;
    _ey[14] = 1.000000;
    _ez[14] = 1.000000;
    _ex[15] = -1.000000;
    _ey[15] = 0.000000;
    _ez[15] = -1.000000;
    _ex[16] = -1.000000;
    _ey[16] = 0.000000;
    _ez[16] = 1.000000;
    _ex[17] = 1.000000;
    _ey[17] = 0.000000;
    _ez[17] = -1.000000;
    _ex[18] = 1.000000;
    _ey[18] = 0.000000;
    _ez[18] = 1.000000;


    _wt[0] = 0.333333;
    _wt[1] = 0.055556;
    _wt[2] = 0.055556;
    _wt[3] = 0.055556;
    _wt[4] = 0.055556;
    _wt[5] = 0.055556;
    _wt[6] = 0.055556;
    _wt[7] = 0.027778;
    _wt[8] = 0.027778;
    _wt[9] = 0.027778;
    _wt[10] = 0.027778;
    _wt[11] = 0.027778;
    _wt[12] = 0.027778;
    _wt[13] = 0.027778;
    _wt[14] = 0.027778;
    _wt[15] = 0.027778;
    _wt[16] = 0.027778;
    _wt[17] = 0.027778;
    _wt[18] = 0.027778;


    {
        // loop over all voxels
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {

                    // natural index for location (i,j)

                    int index = i * ny * nz + j * nz + k;  // column-ordering

                    // initialize density and velocity fields inside the cavity

                    _rho[index] = DENSITY;   // density
                    _ux[index] = 0.0;       // x-component of velocity
                    _uy[index] = 0.0;       // y-component of velocity
                    _uz[index] = 0.0;       // z-component of velocity

                    cpu_ux[index] = 0.0;       // x-component of velocity
                    cpu_uy[index] = 0.0;       // y-component of velocity
                    cpu_uz[index] = 0.0;       // z-component of velocity


                    _sigma[index] = 0.0;       // rate-of-strain field

                    // specify boundary condition for the moving lid

                    //if (j == N - 1) ux[index] = LID_VELOCITY;

                    if ((j == ny - 1)) { _ux[index] = LID_VELOCITY; }

                    // assign initial values for distribution functions
                    // along various aections using equilibriu, functions

                    for (int a = 0; a < Qcc; a++) {

                        int index_f = a + index * Qcc;

                        real edotu = _ex[a] * _ux[index] + _ey[a] * _uy[index] + _ez[a] * _uz[index];
                        real udotu = _ux[index] * _ux[index] + _uy[index] * _uy[index] + _uz[index] * _uz[index];

                        _feq[index_f] = _rho[index] * _wt[a] * (1.0 + 3.0 * edotu + 4.5 * edotu * edotu - 1.5 * udotu);
                       
                        _f[index_f] = _feq[index_f];


                        _f_new[index_f] = _feq[index_f];

                    }

                }
            }
        }
    }


    cudaMemcpy( f, _f, sizeof( real ) * nx * ny * nz * Qcc, cudaMemcpyHostToDevice );
    cudaMemcpy( feq, _feq, sizeof( real ) * nx * ny * nz * Qcc, cudaMemcpyHostToDevice );
    cudaMemcpy( f_new, _f_new, sizeof( real ) * nx * ny * nz * Qcc, cudaMemcpyHostToDevice );

    cudaMemcpy( rho, _rho, sizeof( real ) * nx * ny * nz  , cudaMemcpyHostToDevice );
    cudaMemcpy( ux, _ux, sizeof( real ) * nx * ny * nz  , cudaMemcpyHostToDevice );
    cudaMemcpy( uy, _uy, sizeof( real ) * nx * ny * nz  , cudaMemcpyHostToDevice );
    cudaMemcpy( uz, _uz, sizeof( real ) * nx * ny * nz  , cudaMemcpyHostToDevice );
    cudaMemcpy( sigma, _sigma, sizeof( real ) * nx * ny * nz , cudaMemcpyHostToDevice );


    delete[] _f;
    delete[] _feq;
    delete[] _f_new;
    delete[] _rho;
    delete[] _ux;
    delete[] _uy;
    delete[] _uz;
    delete[] _sigma;
}


void LBMGrid::step( real DENSITY, real  LID_VELOCITY, const real REYNOLDS_NUMBER, bool copy ) {
    this->collideAndStream( DENSITY, LID_VELOCITY, REYNOLDS_NUMBER );
    this->macroVar( DENSITY, LID_VELOCITY, REYNOLDS_NUMBER,copy );
}