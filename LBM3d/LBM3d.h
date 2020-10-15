#pragma once
#pragma once
#include "config3d.h"


void showGraphics( int WIDTH, int HEIGHT, real xmin, real xmax, real ymin, real ymax, const real* ux, const real* uy, const real* uz );


void  D3Q19( real* ex, real* ey, real* ez, int* oppos, real* wt );

void initialize( const int NX, const int NY, const int NZ, const int Q, const real DENSITY, const real LID_VELOCITY,
    real* ex, real* ey, real* ez, int* oppos, real* wt,
    real* rho, real* ux, real* uy, real* uz, real* sigma,
    real* f, real* feq, real* f_new );



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
    real* f_new );      // new distribution function


void macroVar( // READ-ONLY parameters (used by this function but not changed)
    const int NX, const int NY, const int NZ, const int Q, const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER,
    const real* ex, const real* ey, const real* ez, const int* oppos, const real* wt,
    // READ + WRITE parameters (get updated in this function)
    real* rho,         // density
    real* ux,         // X-velocity
    real* uy,         // Y-velocity
    real* uz,         // Y-velocity
    real* sigma,      // rate-of-strain
    real* f,          // distribution function
    real* feq,        // equilibrium distribution function
    real* f_new );    // new distribution function