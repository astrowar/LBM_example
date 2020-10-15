#include "cuda_runtime.h"
#include "lbm_gpu.h"
#include <cstdio>

typedef struct
{
    real x;		//!< x component
    real y;		//!< y component
    real z;		//!< z component
}
vec3_t;

#define Qcc (19)
void   D3Q19( real* gpu_ex, real* gpu_ey, real* gpu_ez, int* gpu_oppos, real* gpu_wt )
{
    real* ex = new real[Qcc];
    real* ey = new real[Qcc];
    real* ez = new real[Qcc];
    int* oppos = new int[Qcc];
    real* wt = new real[Qcc];
    ex[0] = 0.000000;
    ey[0] = 0.000000;
    ez[0] = 0.000000;
    ex[1] = -1.000000;
    ey[1] = 0.000000;
    ez[1] = 0.000000;
    ex[2] = 1.000000;
    ey[2] = 0.000000;
    ez[2] = 0.000000;
    ex[3] = 0.000000;
    ey[3] = -1.000000;
    ez[3] = 0.000000;
    ex[4] = 0.000000;
    ey[4] = 1.000000;
    ez[4] = 0.000000;
    ex[5] = 0.000000;
    ey[5] = 0.000000;
    ez[5] = -1.000000;
    ex[6] = 0.000000;
    ey[6] = 0.000000;
    ez[6] = 1.000000;
    ex[7] = -1.000000;
    ey[7] = -1.000000;
    ez[7] = 0.000000;
    ex[8] = -1.000000;
    ey[8] = 1.000000;
    ez[8] = 0.000000;
    ex[9] = 1.000000;
    ey[9] = -1.000000;
    ez[9] = 0.000000;
    ex[10] = 1.000000;
    ey[10] = 1.000000;
    ez[10] = 0.000000;
    ex[11] = 0.000000;
    ey[11] = -1.000000;
    ez[11] = -1.000000;
    ex[12] = 0.000000;
    ey[12] = -1.000000;
    ez[12] = 1.000000;
    ex[13] = 0.000000;
    ey[13] = 1.000000;
    ez[13] = -1.000000;
    ex[14] = 0.000000;
    ey[14] = 1.000000;
    ez[14] = 1.000000;
    ex[15] = -1.000000;
    ey[15] = 0.000000;
    ez[15] = -1.000000;
    ex[16] = -1.000000;
    ey[16] = 0.000000;
    ez[16] = 1.000000;
    ex[17] = 1.000000;
    ey[17] = 0.000000;
    ez[17] = -1.000000;
    ex[18] = 1.000000;
    ey[18] = 0.000000;
    ez[18] = 1.000000;
    wt[0] = 0.333333;
    wt[1] = 0.055556;
    wt[2] = 0.055556;
    wt[3] = 0.055556;
    wt[4] = 0.055556;
    wt[5] = 0.055556;
    wt[6] = 0.055556;
    wt[7] = 0.027778;
    wt[8] = 0.027778;
    wt[9] = 0.027778;
    wt[10] = 0.027778;
    wt[11] = 0.027778;
    wt[12] = 0.027778;
    wt[13] = 0.027778;
    wt[14] = 0.027778;
    wt[15] = 0.027778;
    wt[16] = 0.027778;
    wt[17] = 0.027778;
    wt[18] = 0.027778;
    oppos[0] = 0;
    oppos[1] = 2;
    oppos[2] = 1;
    oppos[3] = 4;
    oppos[4] = 3;
    oppos[5] = 6;
    oppos[6] = 5;
    oppos[7] = 10;
    oppos[8] = 9;
    oppos[9] = 8;
    oppos[10] = 7;
    oppos[11] = 14;
    oppos[12] = 13;
    oppos[13] = 12;
    oppos[14] = 11;
    oppos[15] = 18;
    oppos[16] = 17;
    oppos[17] = 16;
    oppos[18] = 15;

   // printf( "ex[19] = {" ); for (int i = 0; i < 19; ++i) { printf( "%f ,", ex[i] ); }printf( " };\n" );
   // printf( "ey[19] = {" ); for (int i = 0; i < 19; ++i) { printf( "%f ,", ey[i] ); }printf( " };\n" );
   // printf( "ez[19] = {" ); for (int i = 0; i < 19; ++i) { printf( "%f ,", ez[i] ); }printf( " };\n" );
    //printf( "wt[19] = {" ); for (int i = 0; i < 19; ++i) { printf( "%f ,", wt[i] ); }printf( " };\n" );
    //printf( "oppos[19] = {" ); for (int i = 0; i < 19; ++i) { printf( "%i ,", oppos[i] ); }printf( " };\n" );

  //  /// \brief Velocity vectors for D3Q19 as floating-point values
  //  const vec3_t vel3Dv[19] = {
  //    {  0,  0,  0 },		// zero direction
  //    { -1,  0,  0 },		// 6 directions with velocity 1
  //    {  1,  0,  0 },
  //    {  0, -1,  0 },
  //    {  0,  1,  0 },
  //    {  0,  0, -1 },
  //    {  0,  0,  1 },
  //    { -1, -1,  0 },		// 12 directions with velocity sqrt(2)
  //    { -1,  1,  0 },
  //    {  1, -1,  0 },
  //    {  1,  1,  0 },
  //    {  0, -1, -1 },
  //    {  0, -1,  1 },
  //    {  0,  1, -1 },
  //    {  0,  1,  1 },
  //    { -1,  0, -1 },
  //    { -1,  0,  1 },
  //    {  1,  0, -1 },
  //    {  1,  0,  1 },
  //  };

  //  const real init_weights3D[19] = { (real)1. / 3,
  //(real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18, (real)1. / 18,
  //(real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36, (real)1. / 36 };



  //  for (int i = 0; i < 19; ++i) {
  //      ex[i] = vel3Dv[i].x;
  //      ey[i] = vel3Dv[i].y;
  //      ez[i] = vel3Dv[i].z;
  //      wt[i] = init_weights3D[i];

  //      printf( "ex[%i] = %f ;\n", i, ex[i] );
  //      printf( "ey[%i] = %f ;\n", i, ey[i] );
  //      printf( "ez[%i] = %f ;\n", i, ez[i] );
  //      
  //  }
  //  for (int i = 0; i < 19; ++i)      printf( "wt[%i]  = %f ; \n", i, wt[i] );

  //  // define opposite (anti) aections (useful for implementing bounce back)

  //  const int invVel3D[19] = { 0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15 };
  //  for (int i = 0; i < 19; ++i) {
  //      oppos[i] = invVel3D[i];

  //      printf( "oppos[%i]  = %i ; \n", i, oppos[i] );
  //  }

    cudaMemcpy( gpu_ex, ex, sizeof(real)*Qcc, cudaMemcpyHostToDevice );
    cudaMemcpy( gpu_ey, ey, sizeof( real ) * Qcc, cudaMemcpyHostToDevice );
    cudaMemcpy( gpu_ez, ez, sizeof( real ) * Qcc, cudaMemcpyHostToDevice );
    cudaMemcpy( gpu_oppos, oppos, sizeof( int ) * Qcc, cudaMemcpyHostToDevice );
    cudaMemcpy(  gpu_wt, wt, sizeof( real ) * Qcc, cudaMemcpyHostToDevice );

    delete[] ex;
    delete[] ey;
    delete[] ez;
    delete[] oppos;
    delete[] wt;
}
