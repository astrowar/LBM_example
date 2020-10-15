#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC
#include "config.h"
#include "LBM.h"

#include "..\LBM_gpu\LBM_gpu.h"

int  main( int argc, char* argv[] )
{
    //--------------------------------
    //   Create a WINDOW using GLFW
    //--------------------------------

    LBM_initializeCudaDevice();
    LBM_setup( 256 );

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
    real* f = new real[N * N * Qcc];
    real* feq = new real[N * N * Qcc];
    real* f_new = new real[N * N * Qcc];

    // density and velocity
    real* rho = new real[N * N];
    real* ux = new real[N * N];
    real* uy = new real[N * N];

    // rate-of-strain
    real* sigma = new real[N * N];

    // D3Q9 parameters
    real* ex = new real[Qcc];
    real* ey = new real[Qcc];
    int* oppos = new int[Qcc];
    real* wt = new real[Qcc];

    // fill D3Q9 parameters in constant memory on the GPU
    D2Q9( ex, ey, oppos, wt );

    // launch GPU kernel to initialize all fields
    initialize( N, Qcc, DENSITY, LID_VELOCITY, ex, ey, oppos, wt, rho, ux, uy, sigma, f, feq, f_new );

    // time integration
    int time = 0;
    clock_t t0, tN;
    t0 = clock();

    //---------------------------------------
    // Loop until the user closes the window
    //---------------------------------------

    // specify min and max window coordinates
    real xmin = 0, xmax = N, ymin = 0, ymax = N;


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
        collideAndStream( N, Qcc, DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, ex, ey, oppos, wt, rho, ux, uy, sigma, f, feq, f_new );

        // calculate macroscopic variables
        macroVar( N, Qcc, DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, ex, ey, oppos, wt, rho, ux, uy, sigma, f, feq, f_new );

        // on-the-fly OpenGL graphics
        if (time % 100 == 0)
        {
            showGraphics( WIDTH, HEIGHT, xmin, xmax, ymin, ymax, ux, uy );

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
    delete[] oppos;
    delete[] wt;

    // GLFW clean up
    glfwDestroyWindow( window );
    glfwTerminate();

    // exit main
    return 0;
}
