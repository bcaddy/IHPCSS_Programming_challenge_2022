/**
 * @file main.c
 * @brief This file contains the source code of the application to parallelise.
 * @details This application is a classic heat spread simulation.
 * @author Ludovic Capelli
 **/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <inttypes.h>
#include <math.h>
#include <sched.h>
#include <unistd.h>
#include <string.h>

#include "util.h"

/**
 * @argv[0] Name of the program
 * @argv[1] path to the dataset to load
 **/
int main(int argc, char* argv[])
{
	(void)argc;
	(void)argv;
	MPI_Init(NULL, NULL);

	/////////////////////////////////////////////////////
	// -- PREPARATION 1: COLLECT USEFUL INFORMATION -- //
	/////////////////////////////////////////////////////
	// Ranks for convenience so that we don't throw raw values all over the code
	const int MASTER_PROCESS_RANK = 0;

	// The rank of the MPI process in charge of this instance
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Number of MPI processes in total, commonly called "comm_size" for "communicator size".
	int comm_size;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	/// Rank of the first MPI process
	const int FIRST_PROCESS_RANK = 0;
	/// Rank of the last MPI process
	const int LAST_PROCESS_RANK = comm_size - 1;

	double column_neighbour_rank, row_neighbour_rank;

	// Rank of my up neighbour if any
	switch (my_rank)
	{
	case 0:
		column_neighbour_rank = 1;
		row_neighbour_rank    = 2;
		break;
	case 1:
		column_neighbour_rank = 0;
		row_neighbour_rank    = 3;
		break;
	case 2:
		column_neighbour_rank = 3;
		row_neighbour_rank    = 0;
		break;
	case 3:
		column_neighbour_rank = 2;
		row_neighbour_rank    = 1;
		break;
	}
	//report_placement();

	////////////////////////////////////////////////////////////////////
	// -- PREPARATION 2: INITIALISE TEMPERATURES ON MASTER PROCESS -- //
	////////////////////////////////////////////////////////////////////
	int const rowsPerRank = ROWS/2;
	int const colsPerRank = COLUMNS/2;
	/// Array that will contain my part chunk. It will include the 2 ghost rows (1 up, 1 down)
	double temperatures[rowsPerRank+2][colsPerRank+2];
	/// Temperatures from the previous iteration, same dimensions as the array above.
	double temperatures_last[rowsPerRank+2][colsPerRank+2];
	/// On master process only: contains all temperatures read from input file.
	double all_temperatures[ROWS][COLUMNS];

	// The master MPI process will read a chunk from the file, send it to the corresponding MPI process and repeat until all chunks are read.
	if(my_rank == MASTER_PROCESS_RANK)
	{
		initialise_temperatures(all_temperatures);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	///////////////////////////////////////////
	//     ^                                 //
	//    / \                                //
	//   / | \    CODE FROM HERE IS TIMED    //
	//  /  o  \                              //
	// /_______\                             //
	///////////////////////////////////////////

	////////////////////////////////////////////////////////
	// -- TASK 1: DISTRIBUTE DATA TO ALL MPI PROCESSES -- //
	////////////////////////////////////////////////////////
	double total_time_so_far = 0.0;
	double start_time = MPI_Wtime();

	// Replace with scatter or broadcast
	if(my_rank == MASTER_PROCESS_RANK)
	{
		for(int i = 0; i < comm_size; i++)
		{
			// Is the i'th chunk meant for me, the master MPI process?
			if(i != my_rank)
			{
				// No, so send the corresponding chunk to that MPI process.
				MPI_Ssend(&all_temperatures[i * ROWS_PER_MPI_PROCESS][0], ROWS_PER_MPI_PROCESS * COLUMNS_PER_MPI_PROCESS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
			else
			{
				// Yes, let's copy it straight for the array in which we read the file into.
				// #pragma omp parallel for default(none) shared(all_temperatures, temperatures_last) num_threads(8)
				for(int j = 1; j <= ROWS_PER_MPI_PROCESS; j++)
				{
					for(int k = 0; k < COLUMNS_PER_MPI_PROCESS; k++)
					{
						temperatures_last[j][k] = all_temperatures[j-1][k];
					}
				}
			}
		}
	}
	else
	{
		// Receive my chunk.
		MPI_Recv(&temperatures_last[1][0], ROWS_PER_MPI_PROCESS * COLUMNS_PER_MPI_PROCESS, MPI_DOUBLE, MASTER_PROCESS_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// Copy the temperatures into the current iteration temperature as well
	#pragma omp parallel for default(none) shared(temperatures, temperatures_last)
	for(int i = 1; i <= ROWS_PER_MPI_PROCESS; i++)
	{
		for(int j = 0; j < COLUMNS_PER_MPI_PROCESS; j++)
		{
			temperatures[i][j] = temperatures_last[i][j];
		}
	}

	if(my_rank == MASTER_PROCESS_RANK)
	{
		printf("Data acquisition complete.\n");
	}

	// Wait for everybody to receive their part before we can start processing
	MPI_Barrier(MPI_COMM_WORLD);

	/////////////////////////////
	// TASK 2: DATA PROCESSING //
	/////////////////////////////
	int iteration_count = 0;
	/// Maximum temperature change observed across all MPI processes
	double global_temperature_change;
	/// Maximum temperature change for us
	double my_temperature_change, second_my_temperature_change;
	/// The last snapshot made
	double snapshot[ROWS][COLUMNS];

	MPI_Datatype blockType;
    MPI_Type_vector(rowsPerRank, 1, ROWS, MPI_DOUBLE, &blockType);
    MPI_Type_commit(&blockType);

	MPI_Request requestSendRow;
	MPI_Request requestSendColumn;
	MPI_Request requestRecvRow;
	MPI_Request requestRecvColumn;
	MPI_Request requestGatherSnapshot;

	MPI_Send_init(&temperatures[1][0],                          1,                         blockType , column_neighbour_rank,   0, MPI_COMM_WORLD, &requestUp);
	MPI_Send_init(&temperatures[rowsPerRank][0],       colsPerRank, MPI_DOUBLE, row_neighbour_rank, 0, MPI_COMM_WORLD, &requestDown);

	MPI_Send_init(&temperatures[1][0],   rowsPerRank*colsPerRank, MPI_DOUBLE, MASTER_PROCESS_RANK, 0, MPI_COMM_WORLD, &requestGatherSnapshot);

	MPI_Recv_init(&temperatures_last[rowsPerRank+1][0], 1,                         blockType , column_neighbour_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &requestRecvDown);
	MPI_Recv_init(&temperatures_last[0][0],                      colsPerRank, MPI_DOUBLE, row_neighbour_rank,   MPI_ANY_TAG, MPI_COMM_WORLD, &requestRecvUp);

	while(total_time_so_far < MAX_TIME)
	{
		my_temperature_change = 0.0;

		// ////////////////////////////////////////
		// -- SUBTASK 1: EXCHANGE GHOST CELLS -- //
		// ////////////////////////////////////////

		// Send data to up neighbour for its ghost cells. If my up_neighbour_rank is MPI_PROC_NULL, this MPI_Ssend will do nothing.
		MPI_Start(&requestUp);

		// Send data to down neighbour for its ghost cells. If my down_neighbour_rank is MPI_PROC_NULL, this MPI_Ssend will do nothing.
		MPI_Start(&requestDown);

		/////////////////////////////////////////////
		// -- SUBTASK 2: PROPAGATE TEMPERATURES -- //
		/////////////////////////////////////////////
		my_temperature_change = 0.0; // calculate temperature change

		// Wait for snapshot gather
		if ((iteration_count > 0) && ((iteration_count-1) % SNAPSHOT_INTERVAL == 0)) MPI_Wait(&requestGatherSnapshot, MPI_STATUS_IGNORE);

		#pragma omp parallel for default(none) reduction(max:my_temperature_change) shared(temperatures, temperatures_last)
		for(int i = 2; i <= rowsPerRank-1; i++)
		{
			// Process the cell at the first column, which has no left neighbour
			if(temperatures[i][0] != MAX_TEMPERATURE)
			{
				temperatures[i][0] = (temperatures_last[i-1][0] +
									temperatures_last[i+1][0] +
									temperatures_last[i  ][1]) / 3.0;
				my_temperature_change = fmax(fabs(temperatures[i][0] - temperatures_last[i][0]), my_temperature_change);
			}
			// Process all cells between the first and last columns excluded, which each has both left and right neighbours
			for(int j = 2; j < colsPerRank - 2; j++)
			{
				if(temperatures[i][j] != MAX_TEMPERATURE)
				{
					temperatures[i][j] = 0.25 * (temperatures_last[i-1][j  ] +
												temperatures_last[i+1][j  ] +
												temperatures_last[i  ][j-1] +
												temperatures_last[i  ][j+1]);
					my_temperature_change = fmax(fabs(temperatures[i][j] - temperatures_last[i][j]), my_temperature_change);
				}
			}
			// Process the cell at the last column, which has no right neighbour
			if(temperatures[i][colsPerRank - 1] != MAX_TEMPERATURE)
			{
				temperatures[i][colsPerRank - 1] = (temperatures_last[i-1][colsPerRank - 1] +
																temperatures_last[i+1][colsPerRank - 1] +
																temperatures_last[i  ][colsPerRank - 2]) / 3.0;
				my_temperature_change = fmax(fabs(temperatures[i][colsPerRank - 1]
												- temperatures_last[i][colsPerRank - 1]),
												my_temperature_change);
			}
		}

		// Get the halo data and compute the edge values

		// Receive data from down neighbour to fill our ghost cells. If my down_neighbour_rank is MPI_PROC_NULL, this MPI_Recv will do nothing.
		MPI_Start(&requestRecvDown);

		// Receive data from up neighbour to fill our ghost cells. If my up_neighbour_rank is MPI_PROC_NULL, this MPI_Recv will do nothing.
		MPI_Start(&requestRecvUp);

		// Waits for the MPI sends
		MPI_Wait(&requestUp,   MPI_STATUS_IGNORE);
		MPI_Wait(&requestDown, MPI_STATUS_IGNORE);

		// Waits for the MPI receives
		MPI_Wait(&requestRecvUp,   MPI_STATUS_IGNORE);
		MPI_Wait(&requestRecvDown, MPI_STATUS_IGNORE);

		second_my_temperature_change = 0.0;
		#pragma omp parallel default(none) shared(temperatures, temperatures_last) reduction(max:second_my_temperature_change)
		{
			#pragma omp for
			for(int i = 1; i <= rowsPerRank; i += rowsPerRank-1)
			{
				// Process the cell at the first column, which has no left neighbour
				if(temperatures[i][0] != MAX_TEMPERATURE)
				{
					temperatures[i][0] = (temperatures_last[i-1][0] +
										temperatures_last[i+1][0] +
										temperatures_last[i  ][1]) / 3.0;
					second_my_temperature_change = fmax(fabs(temperatures[i][0] - temperatures_last[i][0]), second_my_temperature_change);
				}
				// Process all cells between the first and last columns excluded, which each has both left and right neighbours
				for(int j = 1; j < colsPerRank - 1; j += colsPerRank - 3)
				{
					if(temperatures[i][j] != MAX_TEMPERATURE)
					{
						temperatures[i][j] = 0.25 * (temperatures_last[i-1][j  ] +
													temperatures_last[i+1][j  ] +
													temperatures_last[i  ][j-1] +
													temperatures_last[i  ][j+1]);
						second_my_temperature_change = fmax(fabs(temperatures[i][j] - temperatures_last[i][j]), second_my_temperature_change);
					}
				}
				// Process the cell at the last column, which has no right neighbour
				if(temperatures[i][colsPerRank - 1] != MAX_TEMPERATURE)
				{
					temperatures[i][colsPerRank - 1] = (temperatures_last[i-1][colsPerRank - 1] +
																	temperatures_last[i+1][colsPerRank - 1] +
																	temperatures_last[i  ][colsPerRank - 2]) / 3.0;
					second_my_temperature_change = fmax(fabs(temperatures[i][colsPerRank - 1]
													- temperatures_last[i][colsPerRank - 1]),
													second_my_temperature_change);
				}
			}

			#pragma omp for
			for(int i = 1; i <= rowsPerRank; i++)
			{
				for(int j = 0; j < colsPerRank; j++)
				{
					temperatures_last[i][j] = temperatures[i][j];
				}
			}
		}

		my_temperature_change = fmax(my_temperature_change, second_my_temperature_change);


		///////////////////////////////////////////////////////
		// -- SUBTASK 3: CALCULATE MAX TEMPERATURE CHANGE -- //
		///////////////////////////////////////////////////////
		// Moved to above in loop

		//////////////////////////////////////////////////////////
		// -- SUBTASK 4: FIND MAX TEMPERATURE CHANGE OVERALL -- //
		//////////////////////////////////////////////////////////
		MPI_Allreduce(&my_temperature_change, &global_temperature_change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		//////////////////////////////////////////////////
		// -- SUBTASK 5: UPDATE LAST ITERATION ARRAY -- //
		//////////////////////////////////////////////////
		// for(int i = 1; i <= ROWS_PER_MPI_PROCESS; i++)
		// {
		// 	for(int j = 0; j < COLUMNS_PER_MPI_PROCESS; j++)
		// 	{
		// 		temperatures_last[i][j] = temperatures[i][j];
		// 	}
		// }

		///////////////////////////////////
		// -- SUBTASK 6: GET SNAPSHOT -- //
		///////////////////////////////////
		if(iteration_count % SNAPSHOT_INTERVAL == 0)
		{
			if(my_rank == MASTER_PROCESS_RANK)
			{
				for(int j = 0; j < comm_size; j++)
				{
					if(j == my_rank)
					{
						// Copy locally my own temperature array in the global one
						#pragma omp parallel for default(none) shared(temperatures, snapshot) firstprivate(j)
						for(int k = 0; k < rowsPerRank; k++)
						{
							for(int l = 0; l < colsPerRank; l++)
							{
								snapshot[j * rowsPerRank + k][l] = temperatures[k + 1][l];
							}
						}
					}
					else
					{
						MPI_Recv(&snapshot[j * rowsPerRank][0], rowsPerRank * colsPerRank, MPI_DOUBLE, j, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
				}

				printf("Iteration %d: %.18f\n", iteration_count, global_temperature_change);
			}
			else
			{
				// Send my array to the master MPI process
				MPI_Start(&requestGatherSnapshot);
			}
		}

		// Calculate the total time spent processing
		if(my_rank == MASTER_PROCESS_RANK)
		{
			total_time_so_far = MPI_Wtime() - start_time;
		}

		// Send total timer to everybody so they too can exit the loop if more than the allowed runtime has elapsed already
		MPI_Bcast(&total_time_so_far, 1, MPI_DOUBLE, MASTER_PROCESS_RANK, MPI_COMM_WORLD);

		// Update the iteration number
		iteration_count++;
	}

	///////////////////////////////////////////////
	//     ^                                     //
	//    / \                                    //
	//   / | \    CODE FROM HERE IS NOT TIMED    //
	//  /  o  \                                  //
	// /_______\                                 //
	///////////////////////////////////////////////

	/////////////////////////////////////////
	// -- FINALISATION 2: PRINT SUMMARY -- //
	/////////////////////////////////////////
	if(my_rank == MASTER_PROCESS_RANK)
	{
		printf("The program took %.2f seconds in total and executed %d iterations.\n", total_time_so_far, iteration_count);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
