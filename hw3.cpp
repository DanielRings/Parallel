#include <iostream> 
#include <stdlib.h>
#include <math.h> 
#include <mpi.h>
#include <string.h>

int main(int argc, char *argv[])
{
   	int rank, size;
   	const double PI 3.14159
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	if(rank == 0)
	{
		if(argc != 3 || size < 1)
		{
			return 1; 
		}
	}

	int N, iterations;
    N = atoi(argv[1]);
	iterations = atoi(argv[2]);

	if(rank == 0)
		std::cout<<"size = "<<size<<"\tN = "<<N<<"\titerations = "<<iterations<<std::endl;

	int nodeRows = N / size;
	int rem = N % size;
	if(rank < rem)
		nodeRows++;
	
	double ** local = new double*[nodeRows];
	for(int i=0; i<nodeRows; i++)
		local[i] = new double[N]; 

	double ** temp = new double*[nodeRows];
	for(int i=0; i<nodeRows; i++)
		temp[i] = new double[N]; 
	

	// initialize
	for(int row = 0; r < nodeRows; row++)
	{
		for(int column = 0; column < N; column++)
		{
			local[row][column] = 0.5;
		}
	}
	if(rank == 0)
	{
		for(int column = 0; column < N; column++)
		{
			local[0][column] = 0;
		}
	}
	if(rank == (size-1))
	{
		for(int column = 0; column < N; column++)
		{
			local[nodeRows-1][column] = 5*sin(PI*pow((((double) column) / N), 2)); 
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double startTime = MPI_Wtime();

	double * lowerSwapRow = new double[N]; 
	double * upperSwapRow = new double[N];

	MPI_Request * requestDown = new MPI_Request;
	MPI_Request * requestUp = new MPI_Request;

	for(int i=0; i < iterations; i++)
	{
		if(rank == 1)
		{	
			std::cout <<"Node: "<<rank<<"\tMADE IT TO STEP: "<<i<<std::endl;
			for(int r=0; r<rowsPerSection; r++){
				for(int c=0; c<N; c++)
					std::cout<<local[r][c]<<" ";
				std::cout<<std::endl;
			}
			std::cout<<std::endl<<std::endl;
		}

		if(size > 1)
		{
			if(rank == 0) // First set of rows
			{
				MPI_Irecv((void*)lowerSwapRow, N, MPI_DOUBLE, rank+1, i, MPI_COMM_WORLD, requestDown);
				MPI_Isend((void*)local[nodeRows-1], N, MPI_DOUBLE, rank+1, i, MPI_COMM_WORLD, requestUp);

				for(int row = 1; row < nodeRows-1; row++)
				{
					temp[row][0] = (local[row-1][N-1]+local[row-1][0]+local[row-1][1]+
								local[row][N-1]+local[row][0]+local[row][1]+
								local[row+1][N-1]+local[row+1][0]+local[row+1][1])/9; 

					for(int column = 1; column < N-1; column++)
					{	
						temp[row][column] = (local[row-1][column-1]+local[row-1][column]+local[row-1][column+1]+
									local[row][column-1]+local[row][column]+local[row][column+1]+
									local[row+1][column-1]+local[row+1][column]+local[row+1][column+1])/9; 
					}	
					
					temp[row][N-1] = (local[row-1][N-2]+local[row-1][N-1]+local[row-1][0]+
								local[row][N-2]+local[row][N-1]+local[row][0]+
								local[row+1][N-2]+local[row+1][N-1]+local[row+1][0])/9; 
				}	

				for(int row = 1; row < nodeRows-2; row++)
				{
					memcpy((void*)local[row], (void*)temp[row], sizeof(double)*N);	
				}

				MPI_Wait(requestDown, MPI_STATUS_IGNORE);
				
				temp[nodeRows-1][0] = (local[nodeRows-2][N-1]+local[nodeRows-2][0]+local[nodeRows-2][1]+
								local[nodeRows-1][N-1]+local[nodeRows-1][0]+local[nodeRows-1][1]+
								lowerSwapRow[N-1]+lowerSwapRow[0]+lowerSwapRow[1])/9; 
				for(int column =1; column < N-1; column++)
				{
					temp[nodeRows-1][column] = (local[nodeRows-2][column-1]+local[nodeRows-2][column]+local[nodeRows-2][column+1]+
								local[nodeRows-1][column-1]+local[nodeRows-1][column]+local[nodeRows-1][column+1]+
								lowerSwapRow[column-1]+lowerSwapRow[column]+lowerSwapRow[column+1])/9; 
				}
				temp[nodeRows-1][N-1] = (local[nodeRows-2][N-2]+local[nodeRows-2][N-1]+local[nodeRows-2][0]+
								local[nodeRows-1][N-2]+local[nodeRows-1][N-1]+local[nodeRows-1][0]+
								lowerSwapRow[N-2]+lowerSwapRow[N-1]+lowerSwapRow[0])/9; 
				
				memcpy((void*)local[nodeRows-2], temp[nodeRows-2], sizeof(double)*N);
				memcpy((void*)local[nodeRows-1], temp[nodeRows-1], sizeof(double)*N);
			}
			else if(rank == (size-1)) // last set of rows
			{
				MPI_Irecv((void*)upperSwapRow, N, MPI_DOUBLE, rank-1, i, MPI_COMM_WORLD, requestUp);
				MPI_Isend((void*)local[0], N, MPI_DOUBLE, rank-1, i, MPI_COMM_WORLD, requestDown);
				
				for(int row = 1; row < nodeRows-1; row++)
				{
					temp[row][0] = (local[row-1][N-1]+local[row-1][0]+local[row-1][1]+
								local[row][N-1]+local[row][0]+local[row][1]+
								local[row+1][N-1]+local[row+1][0]+local[row+1][1])/9; 

					for(int column = 1; column < N-1; column++)
					{	
						temp[row][column] = (local[row-1][column-1]+local[row-1][column]+local[row-1][column+1]+
									local[row][column-1]+local[row][column]+local[row][column+1]+
									local[row+1][column-1]+local[row+1][column]+local[row+1][column+1])/9; 
					}	
					
					temp[row][N-1] = (local[row-1][N-2]+local[row-1][N-1]+local[row-1][0]+
								local[row][N-2]+local[row][N-1]+local[row][0]+
								local[row+1][N-2]+local[row+1][N-1]+local[row+1][0])/9; 

				}	

				for(int row = 2; row < nodeRows-1; row++)
				{
					memcpy((void*)local[row], (void*)temp[row], sizeof(double)*N);	
				}

				MPI_Wait(requestUp, MPI_STATUS_IGNORE);

				temp[0][0] = 		(upperSwapRow[N-1]+upperSwapRow[0]+upperSwapRow[1]+
								local[0][N-1]+local[0][0]+local[0][1]+
								local[1][N-1]+local[1][0]+local[1][1])/9; 
				for(int column =1; column < N-1; column++)
				{
					temp[0][column] = (upperSwapRow[column-1]+upperSwapRow[column]+upperSwapRow[column+1]+
							                     local[0][column-1]+local[0][column]+local[0][column+1]+
							                     local[1][column-1]+local[1][column]+local[1][column+1])/9; 
				}
				temp[0][N-1] =		(upperSwapRow[N-2]+upperSwapRow[N-1]+upperSwapRow[0]+     
						                local[0][N-2]+local[0][N-1]+local[0][0]+
						                local[1][N-2]+local[1][N-1]+local[1][0])/9; 

				memcpy((void*)local[0], temp[0], sizeof(double)*N);
				memcpy((void*)local[1], temp[1], sizeof(double)*N);
			}
			else //every other set of rows
			{
				MPI_Irecv((void*)lowerSwapRow, N, MPI_DOUBLE, rank+1, i, MPI_COMM_WORLD, requestUp);
				MPI_Irecv((void*)upperSwapRow, N, MPI_DOUBLE, rank-1, i, MPI_COMM_WORLD, requestDown);
				MPI_Isend((void*)local[nodeRows-1], N, MPI_DOUBLE, rank+1, i, MPI_COMM_WORLD, requestDown);
				MPI_Isend((void*)local[0], N, MPI_DOUBLE, rank-1, i, MPI_COMM_WORLD, requestUp);
			
				for(int row = 1; row < nodeRows-1; row++)
				{
					temp[row][0] = (local[row-1][N-1]+local[row-1][0]+local[row-1][1]+
								local[row][N-1]+local[row][0]+local[row][1]+
								local[row+1][N-1]+local[row+1][0]+local[row+1][1])/9; 

					for(int column = 1; column < N-1; column++)
					{	
						temp[row][column] = (local[row-1][column-1]+local[row-1][column]+local[row-1][column+1]+
									local[row][column-1]+local[row][column]+local[row][column+1]+
									local[row+1][column-1]+local[row+1][column]+local[row+1][column+1])/9; 
					}	
					
					temp[row][N-1] = (local[row-1][N-2]+local[row-1][N-1]+local[row-1][0]+
								local[row][N-2]+local[row][N-1]+local[row][0]+
								local[row+1][N-2]+local[row+1][N-1]+local[row+1][0])/9; 

				}	

				for(int row = 2; row < nodeRows-2; row++)
				{
					memcpy((void*)local[row], (void*)temp[row], sizeof(double)*N);	
				}

				MPI_Wait(requestDown, MPI_STATUS_IGNORE);
				MPI_Wait(requestUp, MPI_STATUS_IGNORE); 
				
				temp[0][0] = 		(upperSwapRow[N-1]+upperSwapRow[0]+upperSwapRow[1]+
								local[0][N-1]+local[0][0]+local[0][1]+
								local[1][N-1]+local[1][0]+local[1][1])/9; 
				for(int column =1; column < N-1; column++)
				{
					temp[0][column] = (upperSwapRow[column-1]+upperSwapRow[column]+upperSwapRow[column+1]+
							                     local[0][column-1]+local[0][column]+local[0][column+1]+
							                     local[1][column-1]+local[1][column]+local[1][column+1])/9; 
				}
				temp[0][N-1] =		(upperSwapRow[N-2]+upperSwapRow[N-1]+upperSwapRow[0]+     
						                local[0][N-2]+local[0][N-1]+local[0][0]+
						                local[1][N-2]+local[1][N-1]+local[1][0])/9; 

				temp[nodeRows-1][0] = (local[nodeRows-2][N-1]+local[nodeRows-2][0]+local[nodeRows-2][1]+
								local[nodeRows-1][N-1]+local[nodeRows-1][0]+local[nodeRows-1][1]+
								lowerSwapRow[N-1]+lowerSwapRow[0]+lowerSwapRow[1])/9; 

				for(int column =1; column < N-1; column++)
				{
					temp[nodeRows-1][column] = (local[nodeRows-2][column-1]+local[nodeRows-2][column]+local[nodeRows-2][column+1]+
								local[nodeRows-1][column-1]+local[nodeRows-1][column]+local[nodeRows-1][column+1]+
								lowerSwapRow[column-1]+lowerSwapRow[column]+lowerSwapRow[column+1])/9; 
				}
				temp[nodeRows-1][N-1] = (local[nodeRows-2][N-2]+local[nodeRows-2][N-1]+local[nodeRows-2][0]+
								local[nodeRows-1][N-2]+local[nodeRows-1][N-1]+local[nodeRows-1][0]+
								lowerSwapRow[N-2]+lowerSwapRow[N-1]+lowerSwapRow[0])/9; 

				memcpy((void*)local[0], temp[0], sizeof(double)*N);
				memcpy((void*)local[1], temp[1], sizeof(double)*N);
				
				memcpy((void*)local[nodeRows-2], temp[nodeRows-2], sizeof(double)*N);
				memcpy((void*)local[nodeRows-1], temp[nodeRows-1], sizeof(double)*N);
			}
		}

	}


	if(size > 1)
	{
		double localSum = 0; 
		double globalSum = 0; 
		int previousRows = (rank * nodeRows);
		if(rank >= rem){
			previousRows += rem;
		}
		for(int row = 0; row < nodeRows; row++)
		{
			int column = previousRows + row;
			localSum += local[row][column];
		}
		
		MPI_Reduce((void*)&localSum, (void*)&globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	}
	else
	{
		// Fill in lonely case
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double endTime = MPI_Wtime();

	if(rank == 1)
	{
		std::cout << "Local Sum: " << localSum;
		std::cout << "\tGlobal Sum: " << globalSum;
		std::cout << "\tTime: " << (endTime-startTime) << std::endl; 
	}

	MPI_Finalize();    
    	return 0; 
}
    
