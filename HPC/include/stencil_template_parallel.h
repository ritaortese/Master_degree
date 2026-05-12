/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 * See COPYRIGHT in top-level directory.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <math.h>

#include <omp.h>
#include <mpi.h>


#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1

//typedef: create an alias for a data type
typedef unsigned int uint;             //array of integer without sign 

typedef uint    vec2_t[2];             //alias for a two-uint array (it is used to represent 2D vectors)- coordinate (x,y) or dimension (nx,Ny) 
typedef double *restrict buffers_t[4]; //alias for an array of four pointers to double (directions)

typedef struct {   
    double   * restrict data;   //pointer to the data of the plane (with ghost cells)
    vec2_t     size;            //internal dimensions 
} plane_t;                      //structure that represent a plane



// inject_energy(periodic, Nsources, Sources, energy, plane, N)
extern int inject_energy ( const int      ,
                           const int      ,
			               const vec2_t  *,  //array of source coordinates, it's const vec2_t *Sources
			               const double   ,
                                 plane_t *,
                           const vec2_t   );


extern int update_plane ( const int      ,
                          const vec2_t   ,   //size of local patch
                          const plane_t *,   //old grid
                                plane_t * ); //new grid

extern int update_plane_internal(const vec2_t   ,   //size of local patch
                                 const plane_t *,   //old grid
                                       plane_t * ); //update only the internal points of the local patch


extern int update_plane_boundary ( const int      ,
                                     const vec2_t   ,    //size of the grid of MPI tasks
                                     const plane_t *,    //old grid
                                           plane_t * );  //new grid
                                     
extern int get_total_energy( plane_t *,      //pointer to local data and internal dimension
                             double  * );    //total energy 

int initialize ( MPI_Comm *,
                 int       ,
		         int       ,
		         int       ,
		         char    **,
                 vec2_t   *,
                 vec2_t   *,                 
		         int      *,
                 int      *,
		         int      *,
		         int      *,
		         int      *,
		         int      *,
                 vec2_t  **,
                 double   *,
                 plane_t  *,
                 buffers_t *,
                 int *,
                 int * );

// release allocated memory
int memory_release (plane_t   * , buffers_t * );

// print global energy
int output_energy_stat ( int      ,
                         plane_t *,
                         double   ,
                         int      ,
                         MPI_Comm *);


//inject energy at each source location
inline int inject_energy ( const int      periodic,
                           const int      Nsources,
			               const vec2_t  *Sources,  //coordinations of the sources, puntator to an array of couples (x,y)
			               const double   energy,   //energy to inject at each source
                           plane_t *plane,          //local plane where inject energy
                           const vec2_t   N         //grid of MPI tasks, array of 2 integer (number of tasks in x and y direction)
                           )
{
    //calculation of the index in the 1D data array
    register const uint  sizex = plane->size[_x_]+2; //internal dimension of the local patch along x +2 ghost cells
                                                     //take the pointer plane and access the _x_ element of the size array inside the structure


    double * restrict data = plane->data; //pointer to the data of the local plane
    
    #define IDX( i, j ) ( (j)*sizex + (i) ) //macro to calculate the index in the 1D array from 2D coordinates (i,j)
    for (int s = 0; s < Nsources; s++)
        {
            uint x = Sources[s][_x_]; //take the x coordinate of the source s
            uint y = Sources[s][_y_];
            
            data[ IDX(x,y) ] += energy; //inject energy at the source location
            
            if ( periodic )
                {
                    if ( (N[_x_] == 1)  ) //if there is only one MPI task in x direction, if it is periodic there aren't near processes 
                    //the periodicity must be enforced inside the same patch, copying the border values as in the serial version
                    //the sx and dx border are in the same task 
                        {
                            if ( x == 1 )
                                data[IDX(plane->size[_x_]+1, y)] += energy;
                            if ( x == plane->size[_x_] )
                                data[IDX(0, y)] += energy;
                        }
                    
                    if ( (N[_y_] == 1) ) //if there is only one MPI task in y direction, if it is periodic there aren't near processes 
                    //the periodicity must be enforced inside the same patch, copying the border values as in the serial version
                    //the bottom and top border are in the same task 
                        {
                            if ( y == 1 )
                                data[IDX(x, plane->size[_y_]+1)] += energy;
                            if ( y == plane->size[_y_] )
                                data[IDX(x, 0)] += energy;
                        }
                }                
        }
 #undef IDX
    
  return 0;
}





inline int update_plane ( const int      periodic, 
                          const vec2_t   N,         // the grid of MPI tasks
                          const plane_t *oldplane, //plane with the values at time t
                                plane_t *newplane  //plane to write the values at time t+1
                          )
    
{
    register uint  fxsize = oldplane->size[_x_]+2; //width of the local patch +2 ghost cells (total size in x direction)
    //register uint fysize = oldplane->size[_y_]+2; //height of the local patch +2 ghost cells
    
    register uint xsize = oldplane->size[_x_];   //width of the local patch
    register uint  ysize = oldplane->size[_y_];  //height of the local patch
    
    #define IDX( i, j ) ( (j)*fxsize + (i) ) //macro to calculate the index in the 1D array from 2D coordinates (i,j)
    

    double * restrict old = oldplane->data;  //1D array of the old plane
    double * restrict new = newplane->data;  //2D array of the new plane
    
    double alpha = 0.6;
    double constant= 0.25 * (1-alpha);
    #pragma omp parallel for schedule(static) //each thread will update a row block, and for each row it will update all the columns 
    for (uint j = 1; j <= ysize; j++)
        for ( uint i = 1; i <= xsize; i++)
            {

                // NOTE: (i-1,j), (i+1,j), (i,j-1) and (i,j+1) always exist even
                //       if this patch is at some border without periodic conditions;
                //       in that case it is assumed that the +-1 points are outside the
                //       plate and always have a value of 0, i.e. they are an
                //       "infinite sink" of heat

                //explanation of the NOTE: if they are in the patch, they are real cells 
                //otherwise, they are ghost cells that simulate the border conditions
                //with no periodicity, the ghost cells at the border are always 0 (infinite sink of heat)
                
                // five-points stencil formula
                //
                // HINT : check the serial version for some optimization
                //
                

                double result = old[IDX(i,j)] * alpha;

                double sum  = old[IDX(i, j-1)] + old[IDX(i, j+1)] +old[IDX(i-1, j)] + old[IDX(i+1, j)];

                new[IDX(i,j)] = result + sum* constant;
                
            }

    if ( periodic )
        {
            if ( N[_x_] == 1 )
                {
                    for ( uint j = 1; j <= ysize; j++ )
                    {   //first column gets value from right internal column
                        new[IDX(0,j)] = new[IDX(xsize,j)];
                        //right ghost column gets value from left internal column
                        new[IDX(xsize+1,j)] = new[IDX(1,j)];
                    }
                }
  
            if ( N[_y_] == 1 ) 
                {
                    for ( uint i = 1; i <= xsize; i++ )
                    {   //top ghost row gets value from bottom internal row
                        new[IDX(i,0)] = new[IDX(i,ysize)];
                        //bottom ghost row gets value from top internal row
                        new[IDX(i,ysize+1)] = new[IDX(i,1)];
                    }
                    
                }
        }

    
 #undef IDX
  return 0;
}

inline int update_plane_internal(
                          const vec2_t   N,
                          const plane_t *oldplane,
                                plane_t *newplane )
{
    (void)N;

    register uint fxsize = oldplane->size[_x_] + 2; //width of the local patch +2 ghost cells
    register uint xsize  = oldplane->size[_x_];     //width of the local patch
    register uint ysize  = oldplane->size[_y_];     //height of the local patch

    #define IDX(i,j) ((j)*fxsize + (i))

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    double alpha = 0.6;
    double constant= 0.25 * (1-alpha);

    #pragma omp parallel for schedule(static) //each thread will update a row block, and for each row it will update all the columns 
       
    for (uint j = 2; j < ysize; j++)          //update only the internal points
        for (uint i = 2; i < xsize; i++)
        {
            double result = old[IDX(i,j)] * alpha;

            double sum  = old[IDX(i, j-1)] + old[IDX(i, j+1)] +old[IDX(i-1, j)] + old[IDX(i+1, j)];

            new[IDX(i,j)] = result + sum* constant;
        }

    #undef IDX
    return 0;
}

inline int update_plane_boundary(
                          const int      periodic,
                          const vec2_t   N,
                          const plane_t *oldplane,
                                plane_t *newplane )
{
    register uint fxsize = oldplane->size[_x_] + 2;
    register uint xsize  = oldplane->size[_x_];
    register uint ysize  = oldplane->size[_y_];

    #define IDX(i,j) ((j)*fxsize + (i))

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    double alpha = 0.6;
    double constant= 0.25 * (1-alpha);
    uint j_top = 1;
    uint j_bottom = ysize;
    uint i_top = 1;
    uint i_bottom = xsize;

    #pragma omp parallel              //create a parallel region: all the threads will execute the following code, 
                                        
    {
    //top & bottom rows 
    #pragma omp for schedule(static)   //parallelize the loop on the rows, each thread takes a block of rows
    for (uint i = 1; i <= xsize; i++)
    {
        //top row j=1 
        {
            double result = old[IDX(i,j_top)] * alpha;
            
            double sum  = old[IDX(i, j_top-1)] + old[IDX(i, j_top+1)] +old[IDX(i-1, j_top)] + old[IDX(i+1, j_top)];

            new[IDX(i,j_top)] = result + sum* constant;
            
        }

        //bottom row j=ysize 
        {
            double result = old[IDX(i,j_bottom)] * alpha;
            
            double sum  = old[IDX(i, j_bottom-1)] + old[IDX(i, j_bottom+1)] +old[IDX(i-1, j_bottom)] + old[IDX(i+1, j_bottom)];

            new[IDX(i,j_bottom)] = result + sum* constant;
        }
    }


    //left & right columns (no corners) - starting from 2 to avoid the corners that are already updated in the previous loop
    #pragma omp for schedule(static)  //parallelize the loop on the columns, each thread takes a block of columns
    for (uint j = 2; j < ysize; j++)
    {
        //left column 
        {
            double result = old[IDX(i_top,j)] * alpha;
        
            double sum  = old[IDX(i_top, j-1)] + old[IDX(i_top, j+1)] +old[IDX(i_top-1, j)] + old[IDX(i_top+1, j)];

            new[IDX(i_top,j)] = result + sum* constant;
        }

        //right column 
        {
            double result = old[IDX(i_bottom,j)] * alpha;
            
            double sum  = old[IDX(i_bottom, j-1)] + old[IDX(i_bottom, j+1)] +old[IDX(i_bottom-1, j)] + old[IDX(i_bottom+1, j)];

            new[IDX(i_bottom,j)] = result + sum* constant;
        }
    }
}
    /* periodic handling IDENTICO al tuo */
    if (periodic)
    {
        if (N[_x_] == 1)
        {
            #pragma GCC unroll 4  //unroll the loop to improve performance: it will execute 4 iterations of the loop in parallel
            for (uint j = 1; j <= ysize; j++)
            {
                //first column gets value from right internal column
                new[IDX(0,j)] = new[IDX(xsize,j)];
                //right ghost column gets value from left internal column
                new[IDX(xsize+1,j)] = new[IDX(1,j)];
            }
        }

        if (N[_y_] == 1)
        {
            #pragma GCC unroll 4
            for (uint i = 1; i <= xsize; i++)
            {
                //top ghost row gets value from bottom internal row
                new[IDX(i,0)] = new[IDX(i,ysize)];
                //bottom ghost row gets value from top internal row
                new[IDX(i,ysize+1)] = new[IDX(i,1)];
            }
        }
    }

    #undef IDX
    return 0;
}

inline int get_total_energy( plane_t *plane, //local patch
                             double  *energy ) //puntator in which store the total energy
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{

    register const int xsize = plane->size[_x_]; //dimension lecture
    register const int ysize = plane->size[_y_];
    register const int fsize = xsize+2;

    double * restrict data = plane->data; //pointer to the data of the local plane (all)
    
   #define IDX( i, j ) ( (j)*fsize + (i) ) //macro 

   #if defined(LONG_ACCURACY)    
    long double totenergy = 0;
   #else
    double totenergy = 0;    
   #endif

   #pragma omp parallel for reduction(+:totenergy) schedule(static) //parallelize the loop on the rows, and for each row it will sum all the columns, 
                                                                    //the reduction is used to sum the partial results of each thread in a single variable totenergy
    for ( int j = 1; j <= ysize; j++ ) //double loop over the internal cells
        for ( int i = 1; i <= xsize; i++ )
            totenergy += data[ IDX(i, j) ];

    
   #undef IDX

    *energy = (double)totenergy;
    return 0;
}



