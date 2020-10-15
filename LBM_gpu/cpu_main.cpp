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

void showGraphics( int WIDTH, int HEIGHT, real xmin, real xmax, real ymin, real ymax, size_t NX, size_t NY, size_t NZ, const real* ux, const real* uy, const real* uz, real LID_VELOCITY );


class  CursosHandle
{
    double prev_xpos;
    double prev_ypos;
    int state = 0; 
public:
    double rx = 0;
    double ry = 0;

      void cursor_position_callback( GLFWwindow* window  )
    {
          GLdouble xpos, ypos;
          glfwGetCursorPos( window, &xpos, &ypos );
          int w, h;
          glfwGetFramebufferSize( window, &w, &h );
        int mstate = glfwGetMouseButton( window, GLFW_MOUSE_BUTTON_LEFT );
        if (mstate == GLFW_PRESS)
            if (state == 0) {
                prev_xpos = xpos;
                prev_ypos = ypos;
                state = 1;
            }
            else if (state == 1) {
                float dx = (xpos - prev_xpos)/w;
                float dy = (ypos - prev_ypos)/h;

                rx += dx;
                ry += dy;

                prev_xpos = xpos;
                prev_ypos = ypos;

                printf( "%f %f \n", dx, dy );
            }
        if (mstate == GLFW_RELEASE) {
            state = 0;
        }
 
      
    }
};


int  main( int argc, char* argv[] )
{
    //--------------------------------
    //   Create a WINDOW using GLFW
    //--------------------------------

    //LBM_initializeCudaDevice();
   // LBM_setup( 256 );

    const int NX = 32;
    const int NY = 32;
    const int NZ = 64;
    const double  DENSITY = 2.7;            // fluid density in lattice units
    const double  LID_VELOCITY = 0.07;      // lid velocity in lattice units
   // const double  REYNOLDS_NUMBER = 1E6;    // REYNOLDS_NUMBER = LID_VELOCITY * N / kinematicViscosity

    const double  REYNOLDS_NUMBER = 1E6 * (NX/32)  ;    // REYNOLDS_NUMBER = LID_VELOCITY * N / kinematicViscosity


    GLFWwindow* window;
    CursosHandle mouseHandle;
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

    LBMGrid   simulation( NX, NY, NZ );

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
    real xmin = -NX/2, xmax = NX/2, ymin = -NY/2, ymax = NY/2;


    // select background color to be white
    // R = 1, G = 1, B = 1, alpha = 0
    glClearColor( 1.0, 1.0, 1.0, 0.0 );

    // initialize viewing values
    glMatrixMode( GL_PROJECTION );

    // replace current matrix with the identity matrix
    glLoadIdentity();

    // set clipping planes in the X-Y-Z coordinate system
    glOrtho( xmin, xmax, ymin, ymax, -1000.0, 1000.0 );

    glMatrixMode( GL_MODELVIEW );
    //glRotatef( 30.0f, 1, 1, 0 );

 
    
 
     

    while (!glfwWindowShouldClose( window ))
    {
        // increment lattice time
 

        time++;
        simulation.step( DENSITY, LID_VELOCITY, REYNOLDS_NUMBER, time % 10 == 0 );
        // on-the-fly OpenGL graphics
        if (time % 10 == 0)
        {
            showGraphics( WIDTH, HEIGHT, xmin, xmax, ymin, ymax, simulation.nx, simulation.ny, simulation.nz, simulation.cpu_ux, simulation.cpu_uy, simulation.cpu_uz, LID_VELOCITY);

            // swap front and back buffers
            glfwSwapBuffers( window );

            // poll for and processs events
            glfwPollEvents();

            mouseHandle.cursor_position_callback( window );

            glLoadIdentity();
            glRotatef( 30.0f * mouseHandle.ry, 1, 0, 0 );
            glRotatef( 30.0f* mouseHandle.rx, 0, 1, 0 );
            
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

 

    // GLFW clean up
    glfwDestroyWindow( window );
    glfwTerminate();

    // exit main
    return 0;
}
