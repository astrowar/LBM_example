#pragma once

using real = float;

class LBMGrid {
public :

    real* f;
    real* feq;
    real* f_new;

    // density and velocity
    real* rho;
    real* ux;
    real* uy;
    real* uz;

    real* cpu_ux;
    real* cpu_uy;
    real* cpu_uz;

    // rate-of-strain
    real* sigma;

    // D3Q9 parameters
    real* ex;
    real* ey;
    real* ez;
    int* oppos;
    real* wt;

    size_t nx;
    size_t ny;
    size_t nz;

	LBMGrid( size_t nx, size_t ny, size_t nz ); 
    ~LBMGrid();
	void initialize( real DENSITY,real  LID_VELOCITY );
	void step( real DENSITY, real  LID_VELOCITY, const real REYNOLDS_NUMBER,  bool copy );

private:
    void  collideAndStream( const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER );
    void  macroVar( const real DENSITY, const real LID_VELOCITY, const real REYNOLDS_NUMBER, bool copy );
 };

extern "C" {
	int LBM_initializeCudaDevice(); 
	int LBM_setup( unsigned int size );
}
