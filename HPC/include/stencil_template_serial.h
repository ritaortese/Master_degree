/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 * See COPYRIGHT in top-level directory.
 */

//header standard (standard files of the C library)
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <math.h>


//macro definitions: textual substitutions done by the the preprocessor before the compilation
//each time that the preprocessor finds the macro name, it replaces it with the macro definition 
//it see NORTH and it substututes it with 0

//directions
#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

//operations
#define SEND 0
#define RECV 1

//indicate what plane used 
#define OLD 0
#define NEW 1

//coordinates
#define _x_ 0
#define _y_ 1

typedef unsigned int uint;

// ============================================================
//
// function prototypes: functions that the main will call

//PATCH: part of the grid assigned to each MPI task
// initialize the simulation parameters form command line arguments

int initialize ( int      ,
		 char   **,
		 uint    *,
		 int     *,
		 int     *,
		 int     *,
		 int   **,
		 double  *,
		 double **,
                 int     *,
                 int     *,
                 int    * 
		 );

// allocate memory for the simulation planes
int memory_release ( double *, int * );

//extern: make the function visible also from other files
// add energy to the plane at the source locations
extern int inject_energy ( const  int,
                           const int    ,
			   const int   *,
			   const double  ,
			   const uint    [2],
                                 double * );

// update the plane according to the stencil operation  
//using the neighbouring points: 5-point stencil                              
extern int update_plane ( const int       ,
			  const uint    [2],
			  const double   *,
		                double   * );

// compute the total energy in the plane
// to verify the correctness of the simulation and print the total energy of each iteration
extern int get_total_energy( const uint     [2],
                             const double *,
                             double * );


// ============================================================
//
// function definition for inline functions: functions that the compilator can insert directly where they are called
//instead of calling them, to reduce the overhead of function calls

//add energy in the source cells, to simulate points in which heat is injected into the system
//PATCH: part of the grid assigned to each MPI task
inline int inject_energy ( const int     periodic, //1 if periodic boundaries are used
                           const int     Nsources, //number of heat sources
                            const int    *Sources, //array of size 2*Nsources with the x and y coordinates of each source
                            const double  energy,  //amount of energy to inject per source
                            const uint    mysize[2], //size of the local patch (x,y)
                            double *plane ) //array 1D representing the grid
{//macro to convert 2D indices into 1D index(+2 for ghost cells)
    
    register const int fxsize = mysize[_x_]+2;  //it put the value in a register for faster access, than the CPU
    register const int xsize = mysize[_x_];
    register const int ysize = mysize[_y_];


   #define IDX( i, j ) ( (j)*(fxsize) + (i) ) //row j*total width + column i, col e row 0 and mysize+1 are ghost cells
   // value in the grid with ghost cells, linear indexing(not matrix)
    
    for (int s = 0; s < Nsources; s++) {
        
        int x = Sources[2*s]; 
        int y = Sources[2*s+1];
        plane[IDX(x, y)] += energy; //inject energy at source location

        if ( periodic ) //toroidal boundaries 
            {
                if ( x == 1 ) //left column boundary, wrap to right ghost cell
                    plane[IDX(mysize[_x_]+1, y)] += energy; 
                if ( x == xsize ) //right column boundary, wrap to left ghost cell
                    plane[IDX(0, y)] += energy;
                if ( y == 1 ) //bottom row boundary, wrap to bottom ghost cell
                    plane[IDX(x, mysize[_y_]+1)] += energy;
                if ( y == ysize ) //top row boundary, wrap to top ghost cell
                    plane[IDX(x, 0)] += energy;
            }
    }
   #undef IDX //remove MACRO
    
    return 0;
}



// update the plane according to the stencil operation
inline int update_plane ( const int     periodic, 
                          const uint     size[2], //size of the local patch (x,y)
			              const double *old    , //old grid
                          double *new    )       //new grid
/*
 * calculate the new energy values
 * the old plane contains the current data, the new plane
 * will store the updated data
 *
 * NOTE: in parallel, every MPI task will perform the
 *       calculation for its patch
 * Patch: part of the grid assigned to each MPI task
 */
{   //calculation of the util size including ghost cells: ghost: 0, intersnal cells: 1..size, ghost: size+1
    register const int fxsize = size[_x_]+2;  //it put the value in a register for faster access, than the CPU
    //register const int fysize = size[_y_]+2; //we have a square grid so it isn't important
    register const int xsize = size[_x_];
    register const int ysize = size[_y_];
    
    //MACRO: convert 2D indices into 1D index(+2 for ghost cells)
   #define IDX( i, j ) ( (j)*fxsize + (i) )

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    //
    // HINT: in any case, this loop is a good candidate
    //       for openmp parallelization

    double alpha = 0.6; //how much of the heat cell remains in it
    double constant= 0.25 * (1-alpha); //how much of the heat cell is distributed to the neighbors, each neighbor gets 1/4 of the remaining heat
    for (int j = 1; j <= ysize; j++)  //it update every internal cell
        for ( int i = 1; i <= xsize; i++) //there are 2 for so it's sequential
            {
                //
                // five-points stencil formula
                //

                
                // simpler stencil with no explicit diffusivity
                // always conserve the smoohed quantity
                // alpha here mimics how much "easily" the heat
                // travels
                
                double result = old[ IDX(i,j) ] *alpha; //center*alpha: part that remains in the cell
                
                double sum  = old[IDX(i, j-1)] + old[IDX(i, j+1)] +old[IDX(i-1, j)] + old[IDX(i+1, j)];

                new[IDX(i,j)] = result + sum* constant;

                /*   ANOTHER POSSIBLE APPROACH: explicit discretized diffusion

                  // implentation from the derivation of
                  // 3-points 2nd order derivatives (laplacian)
                  // however, that should depends on an adaptive
                  // time-stepping so that given a diffusivity
                  // coefficient the amount of energy diffused is
                  // "small"
                  // however the implicit methods are not stable
                  
               #define alpha_guess 0.5     // mimic the heat diffusivity

                double alpha = alpha_guess;  //now alpha is is a phiysical diffusion coefficient
                double sum = old[IDX(i,j)];  //central value
                
                int   done = 0;
                do  //do.. while: if alpha is too big the simulation can be explode
                    //try with alpha_guess, if the result is acceptable ok otherwise reduce alpha
                    {                
                        double sum_i = alpha * (old[IDX(i-1, j)] + old[IDX(i+1, j)] - 2*sum);
                        double sum_j = alpha * (old[IDX(i, j-1)] + old[IDX(i, j+1)] - 2*sum);
                        result = sum + ( sum_i + sum_j);  //discretization formula
                        double ratio = fabs((result-sum)/(sum!=0? sum : 1.0));  //relative variation/ if sum is =0, divide by 1
                                                                                //fabs: absolute value
                        //ratio: how much the energy in the cell changes
                        done = ( (ratio < 2.0) && (result >= 0) );    // not too fast diffusion and
                                                                     // not so fast that the (i,j)
                                                                     // goes below zero energy (energy becomes negative)
                        alpha /= 2;  //decrease the temporal step
                    }
                while ( !done );
                */

                
            }

    if ( periodic )
        /*
         * propagate boundaries if they are periodic
         *
         * NOTE: when is that needed in distributed memory, if any?
         * 
         * //In MPI this part is needed only in the case inwhich the domain is 1D:
         * N[_x_] = 1, N[_y_] > 1 or N[_x_] > 1, N[_y_] = 1 
         * only 1 row or only 1 column of processes, so the process is neighbor of itself in one direction, 
         * so it needs to update the ghost cells with the values of the internal cells of its patch
         * 
         * In general 2D domain decomposition, ghost cells are updated through 
         * MPI communication with neighboring processes, not by local copying.
         */
        {
            for ( int i = 1; i <= xsize; i++ )
                {
                    new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; //top ghost row gets value from bottom internal row
                    new[ IDX(i, ysize+1) ] = new[ IDX(i, 1) ]; //bottom ghost row gets value from top internal row
                }   
            for ( int j = 1; j <= ysize; j++ )
                {
                    new[ IDX( 0, j) ] = new[ IDX(xsize, j) ]; //right ghost column gets value from left internal column
                    new[ IDX( xsize+1, j) ] = new[ IDX(1, j) ]; //left ghost column gets value from right internal column
                }
        }
    
    return 0;

   #undef IDX
}

 

inline int get_total_energy( const uint     size[2], //local patch
                             const double *plane,   //grid (with ghost cells)
                                   double *energy ) //total energy
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{
    register const int xsize = size[_x_]; //x size
    register const int ysize = size[_y_];
    
    //usual MACRO: convert 2D indices into 1D index(+2 for ghost cells)
   #define IDX( i, j ) ( (j)*(xsize+2) + (i) )

   #if defined(LONG_ACCURACY)    //compilator option: double (faster, less accurate) or long double (slower, more accurate)
                                //if we want more accuracy in the summation
    long double totenergy = 0;
   #else
    double totenergy = 0;    
   #endif

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    for ( int j = 1; j <= ysize; j++ )
        for ( int i = 1; i <= xsize; i++ )
            totenergy += plane[ IDX(i, j) ];  //sum only internal cells, ignoring the ghost cells
    
   #undef IDX

    *energy = (double)totenergy;  //if you used long double, convert to double for output
    return 0;
}
                            
