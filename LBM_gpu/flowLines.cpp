
#include "LBM_gpu.h" 
#include <GLFW/glfw3.h>

#include <math.h>
#include <random>
#include <ctime>
#define NT 40

class ParticleTail {
public:
	float* x;
	float* y;
	float* z;
	int sz;
	ParticleTail() {
		sz = 0;
		x = new float[NT];
		y = new float[NT];
		z = new float[NT];
	}
	void reset() {
		sz = 0;
	}
	void add( float _x, float _y, float _z ) {
		if (sz > NT-1) {
			memcpy( x, x + 1, sz * sizeof( float ) );
			memcpy( y, y + 1, sz * sizeof( float ) );
			memcpy( z, z + 1, sz * sizeof( float ) );
			sz = NT-2;
		}

		x[sz] = _x;
		y[sz] = _y;
		z[sz] = _z;
		sz++;
	}
	
};

class Particles {
public:
	float* x;
	float* y;
	float* z;
	int* dt;
	int n;
	ParticleTail* tails;

	Particles( int _n ) {
		n = _n;
		x = new float[n];
		y = new float[n];
		z = new float[n];

		x = new float[n];
		y = new float[n];
		z = new float[n];

		tails = new ParticleTail[n];

		dt = new int[n];

		init();
		
	}
	void initParticle( int i ) {
		x[i] = 0.001 * (rand() % 1000);
		y[i] = 0.001 * (rand() % 1000);
		z[i] = 0.001 * (rand() % 1000);
		dt[i] = (rand() % 200) + 200;
		tails[i].reset();
	}

	void init(  ) {
		srand( time(NULL) );
		for (int i = 0; i < n; ++i) {
			initParticle( i );
		}
	}

	void update( float du, size_t NX, size_t NY, size_t NZ, const real* ux, const real* uy, const real* uz  ) {

		for (int i = 0; i < n; ++i) {
			int ix = x[i] * NX;
			int iy = y[i] * NY;
			int iz = z[i] * NZ;


			if ((ix < 2 || ix >= NX - 1) || (iy < 2 || iy >= NY - 1) || (iz < 2 || iz >= NZ - 1))
			{
				initParticle( i );
				continue;
			}


			int idx00 = (ix)*NY * NZ + (iy)*NZ + iz;   // point (0,0)

			real vx = ux[idx00] / NX;
			real vy = uy[idx00] / NY;
			real vz = uz[idx00] / NZ;

			x[i] += du * vx;
			y[i] += du * vy;
			z[i] += du * vz;

			tails[i].add( x[i], y[i], z[i] );
			 
			dt[i] = dt[i] - 1;
			 if (dt[i] < 0)initParticle( i );

		}
	}


};


void showLines(   real xmin, real xmax, real ymin, real ymax, real zmin, real zmax,
	size_t NX, size_t NY, size_t NZ,
	const real* ux, const real* uy, const real* uz, real LID_VELOCITY )
{ 
	static Particles *p = new Particles( 200 );
	p->update( 1.0 * NX, NX, NY, NZ, ux, uy, uz );

	glPointSize( 4.0f );
	glColor4f( 0,0,0,1.0f );
	glBegin( GL_POINTS );
	for (int i = 0; i < p->n; ++i) {
		glVertex3f( 
			xmin + (xmax-xmin)* p->x[i],
			ymin + (ymax - ymin) * p->y[i],
			zmin + (zmax - zmin) * (p->z[i] ) );
	}
	glEnd( );

	for (int i = 0; i < p->n; ++i) {
		if (p->tails[i].sz > 1)
		glBegin( GL_LINE_STRIP );
		for (int j = 0; j < p->tails[i].sz ; ++j) {
			glVertex3f( 
				xmin + (xmax - xmin) * 	p->tails[i].x[j], 
				ymin + (ymax - ymin) * 	p->tails[i].y[j], 
				zmin + (zmax - zmin) * 	p->tails[i].z[j] );
		}
		glEnd();
	}


}
