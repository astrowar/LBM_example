// Flow3D.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//

#include <iostream>
#include "lbm_gpu.h"


// LBM3d.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>        // clock_t, clock(), CLOCKS_PER_SEC
 
 

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

    const int NX = 64;
    const int NY = 64;
    const int NZ = 128;
    const double  DENSITY = 2.7;            // fluid density in lattice units
    const double  LID_VELOCITY = 0.07;      // lid velocity in lattice units
    const double  REYNOLDS_NUMBER = 1E6;    // REYNOLDS_NUMBER = LID_VELOCITY * N / kinematicViscosity


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

     LBMGrid   simulation( NX , NY, NZ );
      
     simulation.initialize( DENSITY, LID_VELOCITY );
    // launch GPU kernel to initialize all fields
   
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
        simulation.step( DENSITY, LID_VELOCITY , REYNOLDS_NUMBER );

       
        // on-the-fly OpenGL graphics
        if (time % 20 == 0)
        {
            showGraphics( WIDTH, HEIGHT, xmin, xmax, ymin, ymax, simulation->ux, simulation->uy, simulation->uz );

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
