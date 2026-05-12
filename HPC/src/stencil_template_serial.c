/*
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 */

#define _GNU_SOURCE //activate GNU/POSIX/BSD extensions 

#include "stencil_template_serial.h"
#include <time.h>
#include <sys/stat.h>

//function to write the plane to a binary file for visualization
int dump ( const double *, const uint [2], const char *, double *, double * );

//fucntion to print the grid at each step, only for small grids
void debug_print_grid(double *plane, const uint S[2], int iter)
{
    if (S[_x_] > 10 || S[_y_] > 10)
        return;

    printf("\nGrid after step %d:\n", iter);

    uint fx = S[_x_] + 2;

    for (uint j = 1; j <= S[_y_]; j++) {
        for (uint i = 1; i <= S[_x_]; i++) {
            printf("%6.3f ", plane[j*fx + i]);
        }
        printf("\n");
    }
    printf("\n");
}
// ------------------------------------------------------------------
// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  srand48(0);
  //calculate the wall time(real time) of the simulation, not the CPU time
  struct timespec start, end; 
  clock_gettime(CLOCK_MONOTONIC, &start); //start timer

  int  Niterations;  
  int  periodic;
  uint  S[2];
  int     Nsources;  //number of heat sources
  int    *Sources;  //coordinates of the source energy
  double  energy_per_source; //injected energy per source
  double *planes[2]; //2 pointers to the planes: old and new
  double injected_heat = 0; //total sum of injected heat
  int injection_frequency; //every how many iterations energy is injected 
  int output_energy_at_steps = 0; //if different from 0, print the energy budget at every step
  int debug_grid = 0; //if different from 0, print the grid at every step

  /* argument checking and setting */
  initialize ( argc, argv, &S[0], &periodic, &Niterations,
	       &Nsources, &Sources, &energy_per_source, &planes[0],
	       &output_energy_at_steps, &injection_frequency, &debug_grid );
  
  
  int current = OLD; //active plane 

  if ( injection_frequency > 1 ) //if the injection isn't at every ieration, we need to inject energy at the beginning
    inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current] );
  
  for (int iter = 0; iter < Niterations; iter++) 
    
    {      
      /* new energy from sources */
      
      //if the iteration is multiple of injection_frequency, inject energy
      if ( iter % injection_frequency == 0 ) //every injection_frequency iterations, it adds energy in the source cells
        {
          //modify the current plane
          inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current] );
          injected_heat += Nsources*energy_per_source; //every source adds energy_per_source
        }
                  
      /* update grid points */
      update_plane(periodic, S, planes[current], planes[!current] ); //read from the current plane (old) and write in the opposite one (new)
      //it uses the 5 stencil to update the plane, defined in the header file

      //DEBUG GRID PRINT
      if (debug_grid)
          debug_print_grid(planes[!current], S, iter);

      if ( output_energy_at_steps ) //if different from 0, print the energy budget at every step
        {
          double system_heat; //how much is present form the grid 
          get_total_energy( S, planes[!current], &system_heat); //sum energy in all internal cells
          //it's used to check: system energy vs injected energy
            
          printf("step %d :: injected energy is %g, updated system energy is %g\n", iter, 
          injected_heat, system_heat );

          char filename[100];
          sprintf( filename, "output_serial/plane_%05d.bin", iter );
          
          // Create output directory if it doesn't exist
          if (iter == 0) {
              mkdir("output_serial", 0755);
          }
          
          dump( planes[!current], S, filename, NULL, NULL ); //write the plane to a binary file for visualization
            
        }

      /* swap planes for the new iteration */
      current = !current;  //old is 0, new is 1. the operator ! inverts the binary value
      
    }
  
  
  /* get final heat in the system */
  
  double system_heat;
  get_total_energy( S, planes[current], &system_heat);  //system_heat: total energy in the grid after the last iteration

  printf("injected energy is %g, system energy is %g\n",
	 injected_heat, system_heat );
  
  memory_release( planes[OLD], Sources );  //in reality, planes[OLD] contains both planes
  
  
  clock_gettime(CLOCK_MONOTONIC, &end); //end timer
  double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9; //end-start time in seconds
  printf("Elapsed time: %g seconds\n", elapsed);

  return 0;
}




/* ==========================================================================
   =                                                                        =
   =   initialization                                                       =
   ========================================================================== */



int memory_allocate ( const uint [2],
		      double ** );


int initialize_sources( uint      [2],
			int       ,
			int     **,
      int      );

int initialize ( int      argc,    // the argc (argument count) from command line, how many arguments are passed to the program
		 char   **argv,                // the argv (argument vector) from command line
		 uint     *S,                   // two-uint array defining the x,y dimensions of the grid
		 int     *periodic,            // periodic-boundary tag
		 int     *Niterations,         // how many iterations
		 int     *Nsources,            // how many heat sources
		 int    **Sources,
		 double  *energy_per_source,   // how much heat per source
		 double **planes,
		 int     *output_energy_at_steps,
		 int     *injection_frequency,
     int     *debug_grid
		 )
{
  int ret;
  
  // ··································································
  // set default values

  S[_x_]            = 1000;
  S[_y_]            = 1000;
  *periodic         = 0;
  *Nsources         = 1;
  *Niterations      = 99;
  *output_energy_at_steps = 0;
  *energy_per_source = 1.0;
  *injection_frequency = *Niterations;

  double freq = 0; //relative frequency 
  
  // ··································································
  // process the command line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":x:y:e:E:f:n:p:o:d:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': S[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': S[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'o': *output_energy_at_steps = (atoi(optarg) > 0);
	    break;

	  case 'f': freq = atof(optarg);
	    break;
    
    case 'd': *debug_grid = (atoi(optarg) > 0);
      break;
	    
	  case 'h': printf( "valid options are ( values btw [] are the default values ):\n"
			    "-x    x size of the plate [1000]\n"
			    "-y    y size of the plate [1000]\n"
			    "-e    how many energy sources on the plate [1]\n"
			    "-E    how many energy sources on the plate [1.0]\n"
			    "-f    the frequency of energy injection [0.0]\n"
			    "-n    how many iterations [100]\n"
			    "-p    whether periodic boundaries applies  [0 = false]\n"
			    "-o    whether to print the energy budgest at every step [0 = false]\n"
          "-d    whether to print the grid at every step [0 = false]\n\n"
			    );
	    break;
	    

	  case ':': printf( "option -%c requires an argument\n", optopt);
	    break;

	  case '?': printf(" -------- help unavailable ----------\n");
	    break;
	  }
      }

    if ( opt == -1 )
      break;
  }

  //calculation of injected energy 
  if ( freq == 0 )  //default: inject at every iteration
    *injection_frequency = 1;
  else
    {
      freq = (freq > 1.0 ? 1.0 : freq ); //if freq is greater than 1, set to 1
                                         //otherwise it doesn't change
      *injection_frequency = freq * *Niterations;
    }

  // ··································································
  /*
   * here we should check for all the parms being meaningful
   *
   */

  // 
  /* parameter validation */

  if (S[_x_] <= 0 || S[_y_] <= 0) {
      printf("Error: grid dimensions must be positive\n");
      return 1;
  }

  if (*Niterations <= 0) {
      printf("Error: number of iterations must be > 0\n");
      return 1;
    } else if ( *Niterations > 1000000 ) {
		printf("Error: number of iterations must be less than 1000000\n");
		return 1;
  }

  if (*Nsources <= 0) {
      printf("Error: number of sources must be >= 0\n");
      return 1;
      } else if ((uint) (*Nsources) > S[_x_] * S[_y_] ) {
		printf("Error: number of sources must be less than or equal to the number of grid points\n");
		return 1;
  }

  if (*energy_per_source <= 0) {
      printf("Error: energy per source must be > 0\n");
      return 1;
  }

  if (*injection_frequency <= 0) {
      printf("Error: injection frequency must be > 0\n");
      return 1;
  }

  if ( *output_energy_at_steps != 0 && *output_energy_at_steps != 1 ) {
		printf("Error: output_energy_at_steps must be 0 or 1\n");
		return 1;
	}
  

  // ··································································
  // allocate the needed memory
  //
  ret = memory_allocate( S, planes );
  //planes is a double **, after the call:
  //planes[0] points to the starting address of the block of memory, so to the fisrt plane
  //planes[1] is a offset inside that block, to point to the second plane
  //the block of memory of the 2 planes is the same 
  //so you swap the indices bt old and new, without copying the plane

  // ··································································
  // allocate the heat sources
  //
  ret = initialize_sources( S, *Nsources, Sources, *debug_grid );
  // allocate and initialize Sources with casual positions
  if ( ret != 0 )
    {
      printf("error allocating/initializing the heat sources\n");
      return ret;
    }
  
  return 0;  
}


int memory_allocate ( const uint      size[2],
		            double **planes_ptr )
/*
 * allocate the memory for the planes
 * we need 2 planes: the first contains the
 * current data, the second the updated data
 * in the integration loop then the roles are
 * swapped at every iteration 
 *
 */
{
  if (planes_ptr == NULL )
    // an invalid pointer has been passed
    // manage the situation
    return 1;

  unsigned int bytes = (size[_x_]+2)*(size[_y_]+2); //number of cells of one plane, including ghost cells

  //planes_ptr is an array of 2 pointers to double, so to dynamic memory  
  //sizeof(double): returns the size in bytes of a double variable
  //malloc: ask to the operative system a block of memory of the requested size
  //memset: fill a block of memory with a specific value (0 here)
  
  planes_ptr[OLD] = (double*)malloc( 2*bytes*sizeof(double) ); //allocate a block of memory of 2 times the number of cells of one plane
  memset ( planes_ptr[OLD], 0, 2*bytes*sizeof(double) ); //full the entire block with zeros
  planes_ptr[NEW] = planes_ptr[OLD] + bytes; //it points to the second plane 
      
  return 0;
}


int initialize_sources( uint      size[2],    //dimensions of the grid without ghost cells
                        int       Nsources,   //number of heat sources that you want to create
                        int     **Sources ,   //puntator to the array of coordinates that the function must allocate and fill
                        int debug_grid)       //debug_grid: if different from 0, print the coordinates of the sources
/*
 * randomly spread heat suources
 *
 */
{
  if (debug_grid)  //1 source 
        Nsources = 1; 

  *Sources = (int*)malloc( Nsources * 2 *sizeof(uint) ); //2 integers to store x and y coordinates for each source
  //total elements: Nsources*2
  //it's an int puntator 

  if (debug_grid)
    {
        // 1 source in the middle of the grid, to check if the energy is injected in the right cell
        (*Sources)[0] = size[_x_] / 2;  
        //(*Sources)[0] = 1;             //1 for the periodic test
        (*Sources)[1] = size[_y_] / 2;

    }

  else{ 
  //casual generation for the coordinates : to find the source cell
  for ( int s = 0; s < Nsources; s++ )
    {  
      (*Sources)[s*2] = 1+ lrand48() % size[_x_]; //x coordinate (that is at the 2*s index): from 1 to size[_x_]
      (*Sources)[s*2+1] = 1+ lrand48() % size[_y_];
    }}

  return 0;
}



int memory_release ( double *data, int *sources ) //release allocated memory
  
{
  if( data != NULL )
    free( data );

  if( sources != NULL )
    free( sources );

  
  
  return 0;
}


//write the grid into a binary file for visualization
int dump ( const double *data, const uint size[2], const char *filename, double *min, double *max )
{
  const int xsize = size[0];
  const int ysize = size[1];

  if ( (filename != NULL) && (filename[0] != '\0') ) //file must be valid and not a empty string
    {
      FILE *outfile = fopen( filename, "w" );
      if ( outfile == NULL ){
        printf("error opening file %s for writing\n", filename);
	      return 2;
    }
      fwrite(size, sizeof(uint), 2, outfile); //write the size of the grid at the beginning of the file, to read it later for visualization
      float *array = (float*)malloc( xsize * sizeof(float) ); //temporary array to store one line of the grid
      
      double _min_ = DBL_MAX; //maximum representable double value
      double _max_ = 0; 

      for ( int j = 0; j < ysize; j++ )
	{
	  /*
	  float y = (float)j / size[1];
	  fwrite ( &y, sizeof(float), 1, outfile );
	  */
	  
    //skip first row and column of padding (ghost cells)
	  const double * restrict line = data + (j+1)*(xsize+2) + 1;
	  for ( int i = 0; i < xsize; i++ ) {
	    array[i] = (float)line[i];
	    _min_ = ( line[i] < _min_? line[i] : _min_ ); //update the min and max
	    _max_ = ( line[i] > _max_? line[i] : _max_ ); }
	  
	  fwrite( array, sizeof(float), xsize, outfile );
	}
      
      free( array );
      
      fclose( outfile );

      if ( min != NULL ){   //return min and max of the grid if requested
	      *min = _min_;
      }
      if ( max != NULL ){
	      *max = _max_;
      }
      return 0;
    }

  else return 1;
  
}

