
#include "LBM_gpu.h" 
#include <GLFW/glfw3.h>

#include <math.h>
void showLines( real xmin, real xmax, real ymin, real ymax, real zmin, real zmax,
    size_t NX, size_t NY, size_t NZ,
    const real* ux, const real* uy, const real* uz, real LID_VELOCITY );

void showGraphics( int WIDTH, int HEIGHT,  real xmin, real xmax, real ymin, real ymax, size_t NX, size_t NY, size_t NZ, const real* ux, const real* uy, const real* uz, real LID_VELOCITY )
{
    //--------------------------------
    //  OpenGL initialization stuff 
    //--------------------------------

    //// select background color to be white
    //// R = 1, G = 1, B = 1, alpha = 0
    //glClearColor(1.0, 1.0, 1.0, 0.0);

    //// initialize viewing values
    //glMatrixMode(GL_PROJECTION);

    //// replace current matrix with the identity matrix
    //glLoadIdentity();

    //// set clipping planes in the X-Y-Z coordinate system
    //glOrtho(xmin, xmax, ymin, ymax, -1.0, 1.0);

    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    // clear all pixels
    glClear( GL_COLOR_BUFFER_BIT );

    // 2D array size is identical to the window size in pixels
    const int GNX = WIDTH / 4;
    const int GNY = HEIGHT / 4;

    // calculate pixel size (rectangle to be rendered)
    float dx = (xmax - xmin) / GNX;
    float dy = (ymax - ymin) / GNY;

    float zscl = 1.0;
    if (NX > NY) zscl = float(NZ) / float( NX );
    else  zscl = float( NZ ) / float( NY );

    float zmin = -(xmax - xmin) * zscl/2.0f;
    float zmax = (xmax - xmin) * zscl / 2.0f;
 

    // buffer used to store what we want to plot
    static  float* scalar = new float[WIDTH * HEIGHT];

    // scale factors
    float min_curl = -0.02f;
    float max_curl = 0.02f;

    // loop to fill the buffer that OpenGL will render
    // and assign an appropriate color to that pixel
    glPushMatrix();
    //glTranslatef( 0, 0, -10.0 );
    glDisable( GL_DEPTH_TEST );
    if (zscl > 1.0)
    glScalef(  1.0f/zscl, 1.0f / zscl, 1.0f / zscl );
    showLines(   xmin, xmax, ymin, ymax, zmin, zmax,
        NX, NY, NZ, ux, uy, uz, LID_VELOCITY );

    for (real dz = 0.1; dz <= 0.9; dz = dz + 0.05)
    {
        //glTranslatef( 0, 0, 0.1 );
        for (int i = 0; i < GNX - 1; i++)
        {

            for (int j = 0; j < GNY - 1; j++)
            {

                // map pixel coordinate (i,j) to LBM lattice coordinates (x,y)
                int xin = i * NX / GNX;
                int yin = j * NY / GNY;
                int zin = (1-dz) * NZ;

                if (xin >= NX - 2) continue;
                if (yin >= NY - 2) continue;

                // get locations of 4 data points inside which this pixel lies
                int idx00 = (xin)*NY * NZ + (yin)*NZ + zin;   // point (0,0)
                int idx10 = (xin + 1) * NY * NZ + (yin)*NZ + zin;   // point (1,0)
                int idx01 = (xin)*NY * NZ + (yin + 1) * NZ + zin;   // point (0,1)
                int idx11 = (xin + 1) * NY * NZ + (yin + 1) * NZ + zin;   // point (1,1)

                // additional neighbors for calculating derivatives
                //
                //               0p      1p 
                //               |       |
                //               |       |
                //        m1-----01------11----p1
                //               |       |
                //               |       |
                //               |       |
                //        m0-----00------10----p0
                //               |       |
                //               |       |
                //               0m      1m
                //
                int idxm0 = (xin > 0) ? (xin - 1) * NY * NZ + (yin)*NZ + zin : idx00;
                int idx0m = (yin > 0) ? (xin)*NY * NZ + (yin - 1) * NZ + zin : idx00;
                int idx1m = (yin > 0) ? (xin + 1) * NY * NZ + (yin - 1) * NZ + zin : idx10;
                int idxm1 = (xin > 0) ? (xin - 1) * NY * NZ + (yin + 1) * NZ + zin : idx01;

                int idxp0 = (xin < NX - 2) ? (xin + 2) * NY * NZ + (yin)*NZ + zin : idx10;
                int idxp1 = (xin < NX - 2) ? (xin + 2) * NY * NZ + (yin + 1) * NZ + zin : idx11;
                int idx1p = (yin < NY - 2) ? (xin + 1) * NY * NZ + (yin + 2) * NZ + zin : idx11;
                int idx0p = (yin < NY - 2) ? (xin)*NY * NZ + (yin + 2) * NZ + zin : idx01;


                // calculate the normalized coordinates of the pixel
                float xfl = (float)i * (float)NX / (float)GNX;
                float yfl = (float)j * (float)NY / (float)GNY;
                float x = xfl - (float)xin;
                float y = yfl - (float)yin;

                // calculate "curl" of the velocity field at the 4 data points
                float dVdx_00 = uy[idx10] - uy[idxm0];
                float dVdx_10 = uy[idxp0] - uy[idx00];
                float dVdx_01 = uy[idx11] - uy[idxm1];
                float dVdx_11 = uy[idxp1] - uy[idx01];

                float dUdy_00 = ux[idx01] - ux[idx0m];
                float dUdy_10 = ux[idx11] - ux[idx1m];
                float dUdy_01 = ux[idx0p] - ux[idx00];
                float dUdy_11 = ux[idx1p] - ux[idx10];

                float curl_z_00 = dVdx_00 - dUdy_00;
                float curl_z_10 = dVdx_10 - dUdy_10;
                float curl_z_01 = dVdx_01 - dUdy_01;
                float curl_z_11 = dVdx_11 - dUdy_11;

                // bilinear interpolation
                float ux_interp = ux[idx00] * (1.0 - x) * (1.0 - y) + ux[idx10] * x * (1.0 - y) + ux[idx01] * (1.0 - x) * y + ux[idx11] * x * y;
                float uy_interp = uy[idx00] * (1.0 - x) * (1.0 - y) + uy[idx10] * x * (1.0 - y) + uy[idx01] * (1.0 - x) * y + uy[idx11] * x * y;
                float curl_z_in = curl_z_00 * (1.0 - x) * (1.0 - y) + curl_z_10 * x * (1.0 - y) + curl_z_01 * (1.0 - x) * y + curl_z_11 * x * y;

                // this is the value we want to plot at this pixel (should be in the range [0-1])
                scalar[i * WIDTH + j] = pow( (ux_interp * ux_interp + uy_interp * uy_interp), 0.5f ) / LID_VELOCITY;   // normalized velocity magnitude
               //scalar[i * WIDTH + j] = (max_curl - curl_z_in) / (max_curl - min_curl);                         // normalized vorticity
            }
        }



        for (int i = 0; i < GNX - 1; i++)
        {
            for (int j = 0; j < GNY - 1; j++)
            {

                float x_actual = xmin + i * dx;   // actual x coordinate
                float y_actual = ymin + j * dy;   // actual y coordinate
                float VAL = scalar[i * WIDTH + j];

                float R, G, B, Alpha;
                Alpha = VAL;
                if (VAL < 0.1) continue;
                if (VAL <= 0.5)
                {
                    // yellow to blue transition
                    R = 2 * VAL;
                    G = 2 * VAL;
                    B = 1 - 2 * VAL;
                }
                else
                {
                    // red to yellow transition
                    R = 1;
                    G = 2 - 2 * VAL;
                    B = 0;
                    Alpha = 2 * (VAL - 0.5) + 0.5;
                }

                // rendering the pixel with the appropriate color
                glColor4f( R, G, B,  0.5*Alpha );
               // glRectf( x_actual, y_actual, x_actual + dx, y_actual + dy );

                float z =  zmin + (zmax-zmin)*dz ;
                glBegin( GL_QUADS );
                glVertex3f( x_actual, y_actual, z );
                glVertex3f( x_actual + dx, y_actual, z );

                glVertex3f( x_actual + dx, y_actual + dy, z );
                glVertex3f( x_actual, y_actual + dy, z );

                glEnd();
                
            }
        }
    }
    glPopMatrix();
    // free memory
   // delete[] scalar;
}
