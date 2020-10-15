
#include "cuda_runtime.h"
#include "lbm_gpu.h"

#define Qcc (19)
void   D3Q19( real* ex, real* ey, real* ez, int* oppos, real* wt );

LBMGrid::~LBMGrid()
{
    cudaFree( f );
    cudaFree( feq );
    cudaFree( f_new );

    cudaFree( rho );
    cudaFree( ux );
    cudaFree( uy );
    cudaFree( uz );

    cudaFree( sigma );
    cudaFree( ex );
    cudaFree( ey );
    cudaFree( ez );
    cudaFree( oppos );
    cudaFree( wt );
}
LBMGrid::LBMGrid( size_t _nx, size_t _ny, size_t _nz ) {
    
    this->nx = _nx;
    this->ny = _ny;
    this->nz = _nz;

    cudaMalloc( &f, nx * ny * nz * Qcc * sizeof( real ) );
    cudaMalloc( &feq, nx * ny * nz * Qcc * sizeof( real ) );
    cudaMalloc( &f_new, nx * ny * nz * Qcc * sizeof( real ) );
    
     
   //  f = new real[NX * NY * NZ * Qcc];
   //   feq = new real[NX * NY * NZ * Qcc];
   // f_new = new real[NX * NY * NZ * Qcc];


    cudaMalloc( &rho, nx * ny * nz *  sizeof( real ) );
    cudaMalloc( &ux, nx * ny * nz * sizeof( real ) );
    cudaMalloc( &uy, nx * ny * nz *  sizeof( real ) );
    cudaMalloc( &uz, nx * ny * nz *  sizeof( real ) );

    // density and velocity
    // rho = new real[NX * NY * NZ];
    //  ux = new real[NX * NY * NZ];
   //  uy = new real[NX * NY * NZ];
    //uz = new real[NX * NY * NZ];

     cpu_ux = new real[nx * ny * nz];
     cpu_uy = new real[nx * ny * nz];
     cpu_uz = new real[nx * ny * nz];

    // rate-of-strain
    cudaMalloc( &sigma, nx * ny * nz * sizeof( real ) );
   // sigma = new real[NX * NY * NZ];


    cudaMalloc( &ex, Qcc * sizeof( real ) );
    cudaMalloc( &ey, Qcc * sizeof( real ) );
    cudaMalloc( &ez, Qcc * sizeof( real ) );
    cudaMalloc( &oppos, Qcc * sizeof( int ) );
    cudaMalloc( &wt, Qcc * sizeof( real ) );

    // D3Q9 parameters
    // ex = new real[Qcc];
    //  ey = new real[Qcc];
    //  ez = new real[Qcc];
    //  oppos = new int[Qcc];
    // wt = new real[Qcc];


 
    D3Q19( ex, ey, ez, oppos, wt );

}


 