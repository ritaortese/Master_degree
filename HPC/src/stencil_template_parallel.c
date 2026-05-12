/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */


#include "stencil_template_parallel.h"
#include <float.h>
#include <sys/stat.h>

//function to gather the local planes of all the processes 
double* merge_data(int iter, plane_t *plane, int Rank, int Ntasks, const vec2_t N, const vec2_t S, MPI_Comm *myCOMM_WORLD);

//function to write the plane to a binary file for visualization
int dump ( double *, const uint [2], const char *, double *, double * );



// ------------------------------------------------------------------
// ------------------------------------------------------------------

int main(int argc, char **argv)
{
  MPI_Comm myCOMM_WORLD;            //communicator MPI, crate a local group of tasks to comunicate with different group of them
  int  Rank, Ntasks;                //id of the current MPI task and total number of tasks
  int neighbours[4];               //array to store the neighbours of the current task: 0=N, 1=S, 2=E, 3=W

  int  Niterations;                 //# of iterations
  int  periodic;                    //periodic boundary tag
  vec2_t S, N;                      //size of the global plane, and # of the MPI task in the 2 directions
  
  int      Nsources;                //# of heat sources in the global plane
  int      Nsources_local;          //# of sources in the local patch
  vec2_t  *Sources_local;           //coordinations of the heat sources on the local patch
  double   energy_per_source;       //heat injected at each source

  plane_t   planes[2];              //array of 2 planes (old and new)
  buffers_t buffers[2];             //buffer to communicate with neighbours: buffers[SEND][ direction ], buffers[RECV][ direction ]
  
  int output_energy_stat_perstep;   //if>0 print the statistic 
  int verbose;                      //if>0 print the statistic and dump the plane at every step
  int debug_grid;                    //if>0 print the grid at every step, only for small grids
  /* 

        initialize MPI envrionment 
    
  */
  {
    int level_obtained;
    
    // NOTE: change MPI_FUNNELED if appropriate
    
    //MPI_THREAD_FUNNEDLED: more than 1 thread, but only the master thread can make MPI calls
    //level: the level of thread support provided by the MPI implementation

    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained ); //ask for a level
    if ( level_obtained < MPI_THREAD_FUNNELED ) { //verify of the level of thread support obtained is sufficient
        printf("MPI_thread level obtained is %d instead of %d\n",
        level_obtained, MPI_THREAD_FUNNELED );
        MPI_Finalize();           //close the MPI environment
        exit(1);                  //terminate with an error code
      } 
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);         //assign to Rank the id of the current MPI task
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);       //assign to Ntasks the total number of MPI tasks
    MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD); //create a duplicate of the communicator MPI_COMM_WORLD, so the global one
  } 
  
  double initialization_time = 0; 
  initialization_time -= MPI_Wtime();   //wall-clock time from the starting of the MPI environment, in seconds

  /* 
  
        argument checking and setting 
  
  */
  int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
            neighbours, &Niterations,
            &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
            &planes[0], &buffers[0], &verbose, &debug_grid );
  initialization_time += MPI_Wtime();   //time of initialization: final time - initial time
  
  if ( ret )  //if it is different from 0, it gives an error 
    {
      printf("task %d is opting out with termination code %d\n",
	     Rank, ret );
      
      MPI_Finalize();
      return 0;
    }
  
  
  int current = OLD;
  uint xsize = planes[current].size[_x_];
  uint ysize = planes[current].size[_y_];
  uint xframe = xsize + 2;

 

  double computation_time = 0;
  double communication_time=0;
  double waitall_time=0;
  double total_time = 0;

  total_time -= MPI_Wtime(); //start total timer
  for (int iter = 0; iter < Niterations; ++iter)
    
    { 
      double *data = planes[current].data;
      // Send/Recv: blocked comunication, the task wait until the comunication is completed
      //you can't modify the data in the buffer until the comunication is completed and th eprocess can't do other things 
      
      //Isendv/Irecv: non-blocked comunication, the task can continue its work while the comunication is in progress
      //the function return immediately, without waiting for the comunication to be completed

      MPI_Request reqs[8]; //8 requests: 4 Isend comunication and 4 Irecv comunication 
      int actual_reqs = 0; //counter of the number of requests posted

      /* new energy from sources */
      computation_time -= MPI_Wtime(); 
      inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );
      computation_time += MPI_Wtime(); //time of computation time: final time - start time

      /* -------------------------------------- */

      // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position
      // prepare the bourder data to then update the ghost cells of the neighbour processes

      //if there is a a est neighbour process, fill the SEND buffer with the data of the last column of the local patch
      if (neighbours[EAST] != MPI_PROC_NULL) {
        #pragma GCC unroll 4
        for (uint j = 0; j < ysize; j++) {   //pass all the rows of the local patch
            buffers[SEND][EAST][j] =
                data[(j+1)*xframe + xsize];  // column i = xsize
        }
      }

      //if there is a west neighbour process, fill the SEND buffer with the data of the first column of the local patch
      if (neighbours[WEST] != MPI_PROC_NULL) {
        #pragma GCC unroll 4
        for (uint j = 0; j < ysize; j++) {
            buffers[SEND][WEST][j] =
                data[(j+1)*xframe + 1];   // column i = 1
        }
      }

      //assign a pointer to the correct position of the local patch (without copying)
      //so MPI can directly write the data 
      if (neighbours[NORTH] != MPI_PROC_NULL) {
        buffers[SEND][NORTH] = data + (xframe + 1); //send real row to the last row of the N neighbour
        buffers[RECV][NORTH] = data + 1;            //receive the last row and assign it to the ghost row j=0
    }

      if (neighbours[SOUTH] != MPI_PROC_NULL) {
        buffers[SEND][SOUTH] = data + (ysize * xframe) + 1;
        buffers[RECV][SOUTH] = data + ((ysize+1) * xframe) + 1;
    }



      // [B] perform the halo communications
      //     (1) use Send / Recv
      //    or
      //     (2) use Isend / Irecv
      //         --> can you overlap communication and compution in this way?
      //yes, I can: while the data along the borders are comunicated, I can update the internal points of the patch and then waiting only for update the border cells
      
      communication_time -= MPI_Wtime();

      
      
      #pragma GCC unroll 4
      for (int d = 0; d < 4; d++) {

          if (neighbours[d] != MPI_PROC_NULL){ //if there is neighbour in that direction

          //count the number of elements to send/receive: if the direction is N or S, I need to send/receive xsize elements
          int count = (d < 2 ? (int)xsize : (int)ysize);  // N/S => xsize; E/W => ysize

      
          int recv_tag = d;   //receive from the neighbour 'd' with tag = d
          int send_tag = d^1; //send to the neighbour 'd' with tag = (d^1) that is the opposite direction
                              //d^1 (XOR):swaps 0<->1 and 2<->3

          //non-blocking communication: post the Irecv and Isend
          //the task can continue its work while the comunication is in progress, so the function return immediately
          
          //receive count double from the neighbour 'd' and store it in the RECV buffer
          MPI_Irecv(buffers[RECV][d], count, MPI_DOUBLE,
                    neighbours[d], recv_tag, myCOMM_WORLD, &reqs[actual_reqs++]);

          MPI_Isend(buffers[SEND][d], count, MPI_DOUBLE,
                    neighbours[d], send_tag, myCOMM_WORLD, &reqs[actual_reqs++]);
      
        }
      }

      communication_time += MPI_Wtime(); //time of communication time

      computation_time -= MPI_Wtime(); 
      //while the comunication is in progress, I can update the internal points of the patch
      update_plane_internal(N, &planes[current], &planes[!current]);
      computation_time += MPI_Wtime(); //time of computation time: final time - start time

      waitall_time -= MPI_Wtime();
      //wait until all the non-blocking comunication is completed
      MPI_Waitall(actual_reqs, reqs, MPI_STATUSES_IGNORE);
      waitall_time += MPI_Wtime(); //time of waiting for the comunication to be completed

      
      // [C] copy the haloes data
      //copy the RECV buffers of EAST and WEST to the ghost columns of the local patch, if there are the relative neighbours
      if (neighbours[WEST] != MPI_PROC_NULL) {
        for (uint j = 0; j < ysize; j++) {
            data[(j+1)*xframe + 0] = buffers[RECV][WEST][j];
        }
    }

      if (neighbours[EAST] != MPI_PROC_NULL) {
        for (uint j = 0; j < ysize; j++) {
            data[(j+1)*xframe + (xsize+1)] = buffers[RECV][EAST][j];
        }
    }
      /* --------------------------------------  */
      /* 

            update grid points 
      
      */

      computation_time -= MPI_Wtime();
      //update of the border cells using the data received from the neighbours and stored in the ghost cells
      update_plane_boundary( periodic, N, &planes[current], &planes[!current] );
      //update_plane( periodic, N, &planes[current], &planes[!current] );
      computation_time += MPI_Wtime();

      //DEBUG: PRINT GRID
      if (verbose) {

          double* data = merge_data(iter, &planes[!current], Rank, Ntasks, N, S, &myCOMM_WORLD);

          if (Rank == 0 && data != NULL) {

              printf("\nMerged grid after step %d:\n", iter);

              for (uint j = 0; j < S[_y_]; j++) {
                  for (uint i = 0; i < S[_x_]; i++) {
                      printf("%6.3f ", data[j*S[_x_] + i]);
                  }
                  printf("\n");
              }
              printf("\n");

              free(data);
          }
      }
      /* output if needed */
      if ( output_energy_stat_perstep )
	      output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
      //DEBUG: dump merged global grid
      if (verbose) {
          double *glob = merge_data(iter, &planes[!current], Rank, Ntasks, N, S, &myCOMM_WORLD);
          if (Rank == 0 && glob != NULL) {
              char fname[256];
              // Create output directory if it doesn't exist
              if (iter == 0) {
                  mkdir("output_parallel", 0755);
              }
              sprintf(fname, "output_parallel/parallel_plane_%05d.bin", iter);
              dump(glob, S, fname, NULL, NULL);
              free(glob);
          }
      }
      
      /* swap plane indexes for the new iteration */
      current = !current;
      
    }
  
  total_time += MPI_Wtime();  //total time

  output_energy_stat ( -1, &planes[current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  
  memory_release(planes, buffers);
  
  communication_time += waitall_time; //time of posting the non-blocking comunication + the time of waiting for the comunication to be completed
  double comp_max, waitall_max, comm_max, total_max;
  
  //reduce the times of all the processes and get the maximum time for each phase, so the time of the slowest process: master process receives the results 
  MPI_Reduce(&computation_time, &comp_max, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);   //max computational time among all the processes
  MPI_Reduce(&waitall_time, &waitall_max, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);    //max waiting time among all the processes
  MPI_Reduce(&communication_time, &comm_max, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD); //max communication time among all the processes
  MPI_Reduce(&total_time, &total_max, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);        //max total time among all the processes
  
  if (Rank == 0 || Ntasks == 1) {
    const char *job_name = getenv("JOB_NAME");
    if (!job_name) job_name = "run";

    const char *test_type = getenv("TEST_TYPE");
    if (!test_type) test_type = "undefined";

    const char *nodes = getenv("SLURM_JOB_NUM_NODES");
    if (!nodes) nodes = "1";

    const char *ntasks = getenv("SLURM_NTASKS");
    if (!ntasks) ntasks = "1";

    const char *tasks_per_node = getenv("SLURM_TASKS_PER_NODE");
    if (!tasks_per_node) tasks_per_node = "1";

    const char *cpus_per_task = getenv("SLURM_CPUS_PER_TASK");
    if (!cpus_per_task) cpus_per_task = "1";

    const char *workdir = getenv("SLURM_SUBMIT_DIR");
    if (!workdir) workdir = ".";

    char output_dir[512];
    snprintf(output_dir, sizeof(output_dir), "%s/output_parallel", workdir);

    mkdir(output_dir, 0755);

    char filename[512];
    snprintf(filename, sizeof(filename), "%s/%s.csv", output_dir, job_name);

    int file_exists = (access(filename, F_OK) == 0);

    FILE *f = fopen(filename, "a");
    if (!f) {
        perror("fopen");
        MPI_Abort(myCOMM_WORLD, 1);
    }

    if (!file_exists) {
        fprintf(f,
        "test_type,nodes,ntasks,tasks_per_node,cpus_per_task,Sx,Sy,iterations,total_time,comp_time,comm_time,waitall_time\n");
    }

    fprintf(f,
        "%s,"
        "%s,"
        "%s,"
        "%s,"
        "%s,"
        "%u,%u,"
        "%d,"
        "%.6f,"
        "%.6f,"
        "%.6f,"
        "%.6f\n",
        test_type,
        nodes,
        ntasks,
        tasks_per_node,
        cpus_per_task,
        S[_x_], S[_y_],
        Niterations,
        total_max,
        comp_max,
        comm_max,
        waitall_max
    );

    fclose(f);
}
  MPI_Finalize();
  return 0;
}

/* ==========================================================================
   =                                                                        =
   =   initialization                                                       =
   ========================================================================== */


uint simple_factorization( uint, int *, uint ** );

int initialize_sources( int       ,
			int       ,
			MPI_Comm  *,
			uint      [2],
      vec2_t    ,
      vec2_t    ,
			int       ,
			int      *,
			vec2_t  **, 
      int        );


int memory_allocate ( const int       *,
		            const vec2_t     ,
		            buffers_t *,
		            plane_t   * );
		      

int initialize ( MPI_Comm *Comm,
		 int      Me,                  // the rank of the calling process
		 int      Ntasks,              // the total number of MPI ranks
		 int      argc,                // the argc from command line
		 char   **argv,                // the argv from command line
		 vec2_t  *S,                   // the size of the plane
		 vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
		 int     *periodic,            // periodic-boundary tag
		 int     *output_energy_stat,  // if>0 print the statistic at every step
		 int     *neighbours,          // four-int array that gives back the neighbours of the calling task
		 int     *Niterations,         // how many iterations
		 int     *Nsources,            // how many heat sources
		 int     *Nsources_local,      // how many heat sources in the local patch
		 vec2_t **Sources_local,       // coordinates of the heat sources in the local patch
		 double  *energy_per_source,   // how much heat per source
		 plane_t *planes,              // array of 2 planes (old and new)
		 buffers_t *buffers,            //buffer to communciate
		 int *verbose, 
     int *debug_grid)
{
  int halt = 0; //flag 
  int ret;
  
  
  // ··································································
  // set default values

  (*S)[_x_]         = 10000;
  (*S)[_y_]         = 10000;
  *periodic         = 0;
  *Nsources         = 4;    //number of sources in the global plane 
  *Nsources_local   = 0;    //number of sources in the local patch
  *Sources_local    = NULL; //coordinates of the sources in the local patch, it's a pointer to an array of couples (x,y)
  *Niterations      = 100; //number of iterations
  *energy_per_source = 1.0; //how much energy is injected at each source
  *verbose=0;              //if >0 print the statistic and dump the plane at every step
  *debug_grid=0;           //if >0 print the grid at every step, only for small grids
  *output_energy_stat=0;    //if >0 print the statistic at every step

  if ( planes == NULL ) {  //check if the pointer is valid
    // manage the situation
  }
//set to 0 the sizes of the planes
  planes[OLD].size[0] = planes[OLD].size[1] = 0;  
  planes[NEW].size[0] = planes[NEW].size[1] = 0;
  
  for ( int i = 0; i < 4; i++ )
    neighbours[i] = MPI_PROC_NULL; //initialize the neighbours to null process
    //at the start, you don't know the neighbours of the process
    //then, you will assign the correct values

  //buffers:2x4 matrix- first index the type of buffer (send/recv), second index the direction (N,E,S,W)
  for ( int b = 0; b < 2; b++ )
    for ( int d = 0; d < 4; d++ )
      buffers[b][d] = NULL;  
    //at the start, you don't know what buffers do you need
  
  // ··································································
  // process the command line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:d:")) != -1)
      {
	switch( opt )
	  {
	  case 'x': (*S)[_x_] = (uint)atoi(optarg);
	    break;

	  case 'y': (*S)[_y_] = (uint)atoi(optarg);
	    break;

	  case 'e': *Nsources = atoi(optarg);
	    break;

	  case 'E': *energy_per_source = atof(optarg);
	    break;

	  case 'n': *Niterations = atoi(optarg);
	    break;

	  case 'o': *output_energy_stat = (atoi(optarg) > 0);
	    break;

	  case 'p': *periodic = (atoi(optarg) > 0);
	    break;

	  case 'v': *verbose = atoi(optarg);
	    break;

    case 'd': *debug_grid = (atoi(optarg) > 0);
      break;

	  case 'h': {
	    if ( Me == 0 )
	      printf( "\nvalid options are ( values btw [] are the default values ):\n"
		      "-x    x size of the plate [10000]\n"
		      "-y    y size of the plate [10000]\n"
		      "-e    how many energy sources on the plate [4]\n"
		      "-E    how many energy sources on the plate [1.0]\n"
		      "-n    how many iterations [1000]\n"
		      "-p    whether periodic boundaries applies  [0 = false]\n\n"
      
		      );
	    halt = 1; } //if the help is required, set the halt flag to 1
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

  if ( halt )
    return 1;
  
  
  // ··································································
  /*
   * here we should check that all the parms are meaningful
   *
   */
  if ((*S)[_x_] < 1 || (*S)[_y_] < 1) {
    if (Me == 0)
        printf("Error: grid dimensions must be positive\n");
    return 1;
  }

  if (*Niterations < 1) {
      if (Me == 0)
          printf("Error: number of iterations must be positive\n");
      return 1;
  }

  if (*Nsources < 1) {
      if (Me == 0)
          printf("Error: number of sources must be positive\n");
      return 1;
  } else if ( (uint)(*Nsources) > ((*S)[_x_]*(*S)[_y_]) ) {
      printf("Error: more energy sources than cells in the grid\n");
      return 1;
  }                            

  if (*energy_per_source <= 0.0) {
      if (Me == 0)
          printf("Error: energy per source cannot be negative\n");
      return 1;
  }
  

  
  // ··································································
  /*
   * find a suitable domain decomposition
   * very simple algorithm, you may want to
   * substitute it with a better one
   *
   * the plane Sx x Sy will be solved with a grid
   * of Nx x Ny MPI tasks
   * 
   * grid decomposition in MPI tasks
   */

  int Grid[2]; //grid of MPI tasks
  //formfactor: ratio between the longer and the shorter side of the grid
  //it's similar to 1 id the grid is square, and increases as the grid becomes more rectangular
  double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
  //to decide if the decomposition is 1D or 2D
  int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) ); //Ntasks <= ((int)formfactor+1=1 if you have less processes than the longer side of the grid +1 
  //it is 0 if you have enough processes (2D decomposition ), square

  
  if ( dimensions == 1 ) //if the grid is longer, you divide it in only one axes
    {
      if ( (*S)[_x_] >= (*S)[_y_] )
      	Grid[_x_] = Ntasks, Grid[_y_] = 1;
      else
      	Grid[_x_] = 1, Grid[_y_] = Ntasks;
    }
  else
    {
      int   Nf; //# of factors of Ntasks
      uint *factors; //it finds the factors of Ntasks 
      uint  first = 1; //one of the two dimensions of the grid of processes
      ret = simple_factorization( Ntasks, &Nf, &factors );
      
      for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ )
      	first *= factors[i]; //if the grid is too much rectangular, multiply first by the next factor

      if ( (*S)[_x_] > (*S)[_y_] ) //assign the longer side to the longer side of the grid
      //orient the grid according to the size of the global plane
      	Grid[_x_] = Ntasks/first, Grid[_y_] = first;
      else
      	Grid[_x_] = first, Grid[_y_] = Ntasks/first;
    }

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];
  

  // ··································································
  // my coordinates in the grid of processors
  //
  //Me: rank MPI of the process
  //Grid[_x_]: # of processes in the x direction
  int X = Me % Grid[_x_]; //column of the process
  int Y = Me / Grid[_x_]; //row of the process

  // ··································································
  // find the neighbours of the process in the grid of processes
  //
  //

  if ( Grid[_x_] > 1 )
    {  
      if ( *periodic ) {       
        neighbours[EAST]  = Y*Grid[_x_] + (Me + 1 ) % Grid[_x_]; //(Me + 1 ) % Grid[_x_]: return to 0 when you reach the border
        //Y*Grid[_x_]: shift to the correct row
        //if you are not in the last column, your east neighbour is the next column of the same row
        //otherwise it is the first column of the same row
        neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
        //if you are not in the first column, your west neighbour is Me-1 so the last row

      else {
        neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
        //if you are not in the last column, your east neighbour is the next process Me+1
        //otherwise, it is null process
        neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
        //if you are in the first column, your west neighbour is null process
    }

  if ( Grid[_y_] > 1 )
    {
      if ( *periodic ) {      
          neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
          //if you are in the first row, your north neighbour is the last row
          neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; }
          //if you are in the last row, your south neighbour is the first row
      else {    
          neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
          //if you are not in the first row, your north neighbour is Me - # of processes in x direction
          //otherwise, it is null process
          neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); 
        // if you are not in the last row, your south neighbour is Me + # of processes in x direction
 }
    }

  // ··································································
  // the size of my patch
  //

  /*
   * every MPI task determines the size sx x sy of its own domain
   * REMIND: the computational domain will be embedded into a frame
   *         that is (sx+2) x (sy+2)
   *         the outern frame will be used for halo communication or
   */

  //calculate the size of the local patch
  //divide the global size by the # of processes in that direction
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_]; //internal division for x
  uint r = (*S)[_x_] % Grid[_x_];//remainder
  mysize[_x_] = s + ((uint)X < r); // first r processes get an extra point
  s = (*S)[_y_] / Grid[_y_]; 
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + ((uint)Y < r);

  planes[OLD].size[0] = mysize[0]; //assign the size of the local patch to the plane structure
  //they must be the same
  planes[OLD].size[1] = mysize[1];
  planes[NEW].size[0] = mysize[0];
  planes[NEW].size[1] = mysize[1];
  
//print the decomposition of the grid of processes and the neighbours of each process
  if (*verbose )
    {
      if ( Me == 0 ) { //only the master process prints the grid decomposition
	      printf("Tasks are decomposed in a grid %d x %d\n\n",
		      Grid[_x_], Grid[_y_] );
	      fflush(stdout); //flush the output buffer, so the output is printed immediately
      }

      MPI_Barrier(*Comm); //synchronize all the processes
      
      for ( int t = 0; t < Ntasks; t++ )
      {
        if ( t == Me )//print the info of the current process, in the rank order
          {
            printf("Task %4d :: "
            "\tgrid coordinates : %3d, %3d\n"
            "\tneighbours: N %4d    E %4d    S %4d    W %4d\n",
            Me, X, Y,
            neighbours[NORTH], neighbours[EAST],
            neighbours[SOUTH], neighbours[WEST] );
            fflush(stdout);
          }

        MPI_Barrier(*Comm); //synchronize all the processes 
      }
    }

  
  // ··································································
  // allocate the needed memory
  //
  ret = memory_allocate( neighbours, mysize, buffers, planes );
  

  // ··································································
  // allocate the heat sources
  //
  ret = initialize_sources( Me, Ntasks, Comm, mysize,  *S, *N, *Nsources, Nsources_local, Sources_local, *debug_grid );
  if ( ret!=0 ) {
    printf("task %d is opting out with termination code %d\n", Me, ret);
    memory_release(planes, buffers);
  }
  
  return 0;  
}

//take a # A, find all its factors and store them in factors array
//output: # of factors Nfactors, array of factors of the number A
uint simple_factorization( uint A, int *Nfactors, uint **factors )
/*
 * rought factorization;
 * assumes that A is small, of the order of <~ 10^5 max,
 * since it represents the number of tasks
 #
 */
{
  int N = 0;
  int f = 2;
  uint _A_ = A;

  while ( (uint)f < A )
    {
      while( _A_ % f == 0 ) {  //if f divides _A_ is a factor
        N++;
        _A_ /= f; }

      f++;
    }

  *Nfactors = N;  //now that you know the dimension, you can allocate the array of factors
  uint *_factors_ = (uint*)malloc( N * sizeof(uint) ); //dynamic array of factors

  N   = 0;
  f   = 2;
  _A_ = A;

  while ( (uint)f < A )
    {
      while( _A_ % f == 0 ) {
        _factors_[N++] = f; //fill the array of the factors
        _A_ /= f; }
            f++;
    }

  *factors = _factors_;
  return 0;
}


int initialize_sources( int       Me, //rank MPI of the current process
                        int       Ntasks,
                        MPI_Comm *Comm,
                        vec2_t    mysize,
                        vec2_t S,
                        vec2_t N,
                        int       Nsources,
                        int      *Nsources_local,
                        vec2_t  **Sources,
                        int       debug_grid)

{
  if (debug_grid) {

    //global center
    int global_cx = (int)S[_x_] / 2;
    int global_cy = (int)S[_y_] / 2;

    //coordinates of the grid of processes
    int Grid_x = (int)N[_x_];
    int Grid_y = N[_y_];

    int X = Me % Grid_x;
    int Y = Me / Grid_x;

    //gloal offset of the local patch
    int start_x = 0;
    for (int i = 0; i < X; i++)
        start_x += (int)S[_x_] / Grid_x + (i < ((int)S[_x_] % Grid_x));

    int start_y = 0;
    for (int j = 0; j < Y; j++)
        start_y += (int)S[_y_] / Grid_y + (j < ((int)S[_y_] % Grid_y));

    int local_sx = mysize[_x_];
    int local_sy = mysize[_y_];

    //verify if the global center is in the local patch of the current process
    if (global_cx >= start_x &&
        global_cx < start_x + local_sx &&
        global_cy >= start_y &&
        global_cy < start_y + local_sy) {

        *Nsources_local = 1;

        vec2_t *helper = malloc(sizeof(vec2_t));

        //store the coordinates of the global center in the local coordinates of the patch
        helper[0][_x_] = global_cx - start_x + 1;
        helper[0][_y_] = global_cy - start_y + 1;

        *Sources = helper;
    }
    else {
        *Nsources_local = 0;
        *Sources = NULL;
    }

    return 0;
}
  srand48(time(NULL) ^ Me); //random seed for each process
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  
  if ( Me == 0 )
    {
      for ( int i = 0; i < Nsources; i++ )
        tasks_with_sources[i] = (int)lrand48() % Ntasks;
        //assign randomly the sources to the processes, 
        //so you have an array of size Nsources, where each element is the rank of the process that has a source
    }
  
  MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm ); //broadcast the array of sources to all the processes

  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ )
    nlocal += (tasks_with_sources[i] == Me); //count how many sources are assigned to the current process
  *Nsources_local = nlocal;
  
  if ( nlocal > 0 ) //if the process has sources, it generates local coordinates for the sources
    {
      vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
      for ( int s = 0; s < nlocal; s++ )
        {
          helper[s][_x_] = 1 + lrand48() % mysize[_x_];
          helper[s][_y_] = 1 + lrand48() % mysize[_y_];
        }

      *Sources = helper;
    }
  
  free( tasks_with_sources );

  return 0;
}



int memory_allocate ( 
                      const int       *neighbours  , //neighbours of task 
		                  const vec2_t     N           , //grid size of MPI tasks
                      buffers_t *buffers_ptr ,       //pointers to the buffers for communication
                      plane_t   *planes_ptr          //pointer to the planes for computation
		      )

{ (void)N;
    /*
      here you allocate the memory buffers that you need to
      (i)  hold the results of your computation
      (ii) communicate with your neighbours

     */

  if (planes_ptr == NULL )
    {
      printf("Invalid pointer: planes_ptr is NULL\n");
      return 1;
    }
    


  if (buffers_ptr == NULL )
    {
      printf("Invalid pointer: buffers_ptr is NULL\n");
      return 1;
    }
    
    

  // ··················································
  // allocate memory for data 
  
  // Each MPI task owns a local patch of the global grid.
  //
  // For each task we allocate two planes:
  //   - OLD : contains the solution at time t
  //   - NEW : will contain the updated solution at time t+1
  //
  // Each plane includes an additional one-cell-wide halo
  // (ghost layer) on all sides. These ghost cells are used
  // to store data received from neighbouring MPI tasks
  // during halo exchange.
  unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2);

  // allocate the memory for the "OLD" plane
  planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) ); //ask to OS a block of memory
  if ( planes_ptr[OLD].data == NULL )
    {
      printf("Malloc fail: planes_ptr[OLD].data is NULL\n");
      return 1;
    }
  memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) ); //full the block of memory with zeros

  // allocate the memory for the "NEW" plane
  planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
  if ( planes_ptr[NEW].data == NULL )
    {
      printf("Malloc fail: planes_ptr[NEW].data is NULL\n");
      return 1;
    }
  memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );


  // ··················································
  // allocate buffers
  //
  
  // We allocate the buffers for EAST and WEST communication, that are not contiguous in memory, 
  // to send and receive the columns of the local patch to/from the neighbours.

  // For NORTH and SOUTH no separate buffers are required,
  // since rows are already contiguous in memory and can be sent directly.
  
  // We have buffers[2][4]:
  //  - first index: 0 for send, 1 for receive
  //  - second index: direction of the communication (NORTH, EAST, SOUTH, WEST)
 
  uint ysize = planes_ptr[OLD].size[_y_];

  //EAST buffers
  if (neighbours[EAST] != MPI_PROC_NULL) { //if there is an east neighbour process, you allocate the east buffers

    buffers_ptr[SEND][EAST] = malloc(ysize * sizeof(double));
    if (buffers_ptr[SEND][EAST] == NULL) {
        printf("Malloc failed for SEND EAST buffer\n");
        return 1;
    }

    buffers_ptr[RECV][EAST] = malloc(ysize * sizeof(double));
    if (buffers_ptr[RECV][EAST] == NULL) {
        printf("Malloc failed for RECV EAST buffer\n");
        return 1;
    }

  } else {
      buffers_ptr[SEND][EAST] = NULL;
      buffers_ptr[RECV][EAST] = NULL;
  }

  //WEST buffers
  if (neighbours[WEST] != MPI_PROC_NULL) { //if there is a west neighbour process, you allocate the west buffers

    buffers_ptr[SEND][WEST] = malloc(ysize * sizeof(double));
    if (buffers_ptr[SEND][WEST] == NULL) {
        printf("Malloc failed for SEND WEST buffer\n");
        return 1;
    }

    buffers_ptr[RECV][WEST] = malloc(ysize * sizeof(double));
    if (buffers_ptr[RECV][WEST] == NULL) {
        printf("Malloc failed for RECV WEST buffer\n");
        return 1;
    }

  } else {
      buffers_ptr[SEND][WEST] = NULL;
      buffers_ptr[RECV][WEST] = NULL;
  }

  //NORTH and SOUTH: NO allocation required
  buffers_ptr[SEND][NORTH] = NULL; 
  buffers_ptr[RECV][NORTH] = NULL; 

  buffers_ptr[SEND][SOUTH] = NULL; 
  buffers_ptr[RECV][SOUTH] = NULL;

  // ··················································
  
  
  return 0;
}



 int memory_release ( plane_t   *planes,
                      buffers_t *buffers
		     )
  
{
  //free the memory allocated for the planes and the buffers

  if ( planes != NULL )
    {
      if ( planes[OLD].data != NULL )
	      free (planes[OLD].data);
      
      if ( planes[NEW].data != NULL )
	      free (planes[NEW].data);
    }
  
  if ( buffers != NULL )
    {
      if ( buffers[SEND][EAST] != NULL )
        free (buffers[SEND][EAST]);

      if ( buffers[RECV][EAST] != NULL )
        free (buffers[RECV][EAST]);

      if ( buffers[SEND][WEST] != NULL )
        free (buffers[SEND][WEST]);

      if ( buffers[RECV][WEST] != NULL )
        free (buffers[RECV][WEST]);

    }
      
  return 0;
}



int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm )
{

  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  
  if ( Me == 0 )
    {
      if ( step >= 0 ){
        printf(" [ step %4d ] ", step ); 
        fflush(stdout);
        }

      
      printf( "total injected energy is %g, "
	      "system energy is %g "
	      "( in avg %g per grid point)\n",
	      budget,
	      tot_system_energy,
	      tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
    }
  
  return 0;
}

double* merge_data(int iter, 
                   plane_t *plane,
                   int Rank,
                   int Ntasks,
                   const vec2_t N,   // grid MPI
                   const vec2_t S,   //global dimension
                   MPI_Comm *myCOMM_WORLD)
{
    (void)iter;
  //reconstruct the global grid (without ghost cells) on the master process from the local patchs of each process,
  //so it can be writeed on a file for visualization
    uint sizeX = plane->size[_x_];
    uint sizeY = plane->size[_y_];
    uint frameX = sizeX + 2;

    size_t local_size = sizeX * sizeY;

    //copy the local patch (without ghost cells) into a contiguous buffer 
    double *localbuf = malloc(local_size * sizeof(double));
    for (uint j = 0; j < sizeY; j++)
        for (uint i = 0; i < sizeX; i++)
            localbuf[j*sizeX + i] =
                plane->data[(j+1)*frameX + (i+1)];

    if (Ntasks == 1)
        return localbuf;

    int Grid_x = N[_x_];   //number of processes in x direction
    int Grid_y = N[_y_];   //number of processes in y direction

    int global_x = S[_x_]; //size of the global grid
    int global_y = S[_y_]; //size of the global grid

    double *gathered = NULL; //buffer to gather the local patches on the master process
    int *recvcounts = NULL;  //count of the real cells sent by each process
    int *displs = NULL;      //displacements of the local patches in the global grid

    recvcounts = (int*)calloc((size_t)Ntasks, sizeof(int));
    displs = (int*)calloc((size_t)Ntasks, sizeof(int));
    if (Rank == 0) {
        gathered = (double*)malloc(global_x*global_y*sizeof(double));

        int offset = 0;
        for (int r = 0; r < Ntasks; r++) {
            int gx = r % Grid_x;
            int gy = r / Grid_x;

            int sx = global_x / Grid_x + (gx < (global_x % Grid_x));
            int sy = global_y / Grid_y + (gy < (global_y % Grid_y));

            recvcounts[r] = sx * sy;
            displs[r] = offset;
            offset += recvcounts[r];
        }
    }

    //gather the local patches on the master process, it contains:
    //all the blocks of the ranks, ordered by rank, but not in the correct geometric position
    MPI_Gatherv(localbuf, local_size, MPI_DOUBLE,
                gathered, recvcounts, displs,
                MPI_DOUBLE, 0, *myCOMM_WORLD);

    free(localbuf);

    if (Rank != 0)
        return NULL;

    ///reconstruction of the correct geometric position, so the global grid
    double *global_grid = malloc(global_x*global_y*sizeof(double));

    for (int r = 0; r < Ntasks; r++) {  //for each rank 
        
        //coordinates of the block 
        int gx = r % Grid_x; 
        int gy = r / Grid_x;

        //size of the patch
        int sx = global_x / Grid_x + (gx < (global_x % Grid_x));
        int sy = global_y / Grid_y + (gy < (global_y % Grid_y));

        //position of the block 
        int px = 0; //offset x
        for (int i = 0; i < gx; i++)
            px += global_x / Grid_x + (i < (global_x % Grid_x));

        int py = 0; //offset y
        for (int j = 0; j < gy; j++)
            py += global_y / Grid_y + (j < (global_y % Grid_y));

        //copy the block sx x sy from the gathered buffer to the correct position in the global grid
        int base = displs[r];

        for (int j = 0; j < sy; j++)
            for (int i = 0; i < sx; i++)
                global_grid[(py+j)*global_x + (px+i)] =
                    gathered[base + j*sx + i];
    }

    free(gathered);
    free(recvcounts);
    free(displs);

    return global_grid;
}


//write the grid into a binary file for visualization
int dump ( double *data, const uint size[2], const char *filename, double *min, double *max )
{
    uint xsize = size[0];
    uint ysize = size[1];
    //uint xframe = xsize + 2; //size of the local patch + 2 for the halo
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

      for ( uint j = 0; j < ysize; j++ )
	{
	  /*
	  float y = (float)j / size[1];
	  fwrite ( &y, sizeof(float), 1, outfile );
	  */
	  
    //skip first row and column of padding (ghost cells)
	  const double * restrict line = data + j*xsize;
	  for ( uint i = 0; i < xsize; i++ ) {
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

