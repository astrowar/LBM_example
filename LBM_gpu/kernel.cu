
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "LBM_gpu.h"

#include <stdio.h>
#include <iostream>
#include <math.h>
 



int LBM_initializeCudaDevice()
{
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice( 0 );
    if (cudaStatus != cudaSuccess) {
        fprintf( stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?" );
        goto Error;
    }

Error:

    return cudaStatus;

}

int LBM_setup( unsigned int size )
{
    cudaError_t cudaStatus;
    int* dev_a = 0;
    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc( (void**)&dev_a, size * sizeof( int ) );
    if (cudaStatus != cudaSuccess) {
        fprintf( stderr, "cudaMalloc failed!" );
        goto Error;
    }


Error:
    return cudaStatus;
}

 



// Kernel function to add the elements of two arrays
__global__
void add( int n, float* x, float* y )
{
    for (int i = 0; i < n; i++)
        y[i] = x[i] + y[i];
}

int not_main( void )
{
    int N = 1024;
    float* x, * y;

    // Allocate Unified Memory – accessible from CPU or GPU
    cudaMallocManaged( &x, N * sizeof( float ) );
    cudaMallocManaged( &y, N * sizeof( float ) );

    // initialize x and y arrays on the host
    for (int i = 0; i < N; i++) {
        x[i] = 1.0f;
        y[i] = 2.0f;
    }

    // Run kernel on 1M elements on the GPU
    add << <1, 1 >> > (N, x, y);

    // Wait for GPU to finish before accessing on host
    cudaDeviceSynchronize();

    // Check for errors (all values should be 3.0f)
    float maxError = 0.0f;
    for (int i = 0; i < N; i++)
        maxError = fmax( maxError, fabs( y[i] - 3.0f ) );
    std::cout << "Max error: " << maxError << std::endl;

    // Free memory
    cudaFree( x );
    cudaFree( y );

    return 0;
}

 