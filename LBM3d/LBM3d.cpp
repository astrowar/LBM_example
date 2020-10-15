// LBM3d.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC
#include "config3d.h"
#include "LBM3d.h"

extern "C" {
    int LBM_initializeCudaDevice();
    int LBM_setup( unsigned int size );
}

int  main( int argc, char* argv[] )
{
    //--------------------------------
    //   Create a WINDOW using GLFW
    //--------------------------------

    //LBM_initializeCudaDevice();
   // LBM_setup( 256 );

    GLFWwindow* window;

    // initialize the library
    if (!glfwInit())
        return -1;

    // window size for displaying graphics
    int WIDTH = 800;
    int HEIGHT = 800;

    // set the window's display mode
    window = glfwCreateWindow( WIDTH, HEIGHT, "Flow inside a square cavity", NULL, NULL );
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // make the windows context current
    glfwMakeContextCurrent( window );

    // allocate memory

    // distribution functions
    real* f = new real[NX *NY * NZ * Qcc];
    real* feq = new real[NX * NY * NZ * Qcc];
    real* f_new = new real[NX * NY * NZ * Qcc];

    // density and velocity
    real* rho = new real[NX * NY * NZ];
    real* ux = new real[NX * NY * NZ];
    real* uy = new real[NX * NY * NZ];
    real* uz = new real[NX * NY * NZ];

    // rate-of-strain
    real* sigma = new real[NX * NY * NZ];

    // D3Q9 parameters
    real* ex = new real[Qcc];
    real* ey = new real[Qcc];
    real* ez = new real[Qcc];
    int* oppos = new int[Qcc];
    real* wt = new real[Qcc];

    // fill D3Q9 parameters in constant memory on the GPU
    D3Q19( ex, ey, ez, oppos, wt );

    // launch GPU kernel to initialize all fields
    initialize( NX,NY,NZ, Qcc, DENSITY, LID_VELOCITY, ex, ey, ez,oppos, wt, rho, ux, uy, uz, sigma, f, feq, f_new );

    // time integration
    int time = 0;
    clock_t t0, tN;
    t0 = clock();

    //---------------------------------------
    // Loop until the user closes the window
    //---------------------------------------

    // specify min and max window coordinates
    real xmin = 0, xmax = NX, ymin = 0, ymax = NY;


    // select background color to be white
    // R = 1, G = 1, B = 1, alpha = 0
    glClearColor( 1.0, 1.0, 1.0, 0.0 );

    // initialize viewing values
    glMatrixMode( GL_PROJECTION );

    // replace current matrix with the identity matrix
    glLoadIdentity();

    // set clipping planes in the X-Y-Z coordinate system
    glOrtho( xmin, xmax, ymin, ymax, -1.0, 1.0 );



    while (!glfwWindowShouldClose( window ))
    {
        // increment lattice time
        time++;

        // collision and streaming
        collideAndStream( NX, NY, NZ,  Qcc, DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, ex, ey, ez, oppos, wt, rho, ux, uy,uz, sigma, f, feq, f_new );

        // calculate macroscopic variables
        macroVar( NX, NY, NZ, Qcc, DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, ex, ey,ez, oppos, wt, rho, ux, uy, uz, sigma, f, feq, f_new );

        // on-the-fly OpenGL graphics
        if (time % 20 == 0)
        {
            showGraphics( WIDTH, HEIGHT, xmin, xmax, ymin, ymax, ux, uy,uz );

            // swap front and back buffers
            glfwSwapBuffers( window );

            // poll for and processs events
            glfwPollEvents();
        }

        // calculate and print the number of lattice time-steps per second
        tN = clock() - t0;
        if (0) {
            std::cout << "Lattice time " << time
                << " clock ticks " << tN
                << " wall clock time " << tN / CLOCKS_PER_SEC
                << " lattice time steps per second = " << (float)CLOCKS_PER_SEC * time / (float)tN
                << std::endl;
        }
    }

    // free memory for LBM buffers
    delete[] f;
    delete[] feq;
    delete[] f_new;
    delete[] rho;
    delete[] ux;
    delete[] uy;
    delete[] sigma;
    delete[] ex;
    delete[] ey;
    delete[] ez;
    delete[] oppos;
    delete[] wt;

    // GLFW clean up
    glfwDestroyWindow( window );
    glfwTerminate();

    // exit main
    return 0;
}
