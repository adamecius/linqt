
#ifndef PARALLELIZATION_SCHEME_HPP
#define PARALLELIZATION_SCHEME_HPP

#include <fstream>
#include <sstream>
#include <iostream>
#include "types_definitions.hpp"
#include "string_util.h"
#include <mpi.h>
#include "mpi_ofile.hpp"
#include <sstream>
#include "parser.hpp"
#include <cstring>
#include <sys/time.h>
#include <ctime>
#include <cassert>


void GatherAndReduce(	kpm::real* recv, int RecvSize,
						kpm::real* send, int SendSize,
						int color, int key,
						const MPI_Comm &outer_comm)
		{
			kpm::real* trecv= NULL;
			int row_root=0,row_rank;
			//Divide the total list of processes into rows
			MPI_Comm  row_comm;
			MPI_Comm_split(outer_comm, color , key , &row_comm);
			MPI_Comm_rank(row_comm,&row_rank);
			//Initialize the receiver only at the root of each color
			if ( row_rank== row_root )
				trecv= new kpm::real[RecvSize];


			//Gather the data only at the root  
			MPI_Gather(	&send[0] ,SendSize, MPI_DOUBLE, 
						&trecv[key],SendSize, MPI_DOUBLE,
						row_root , row_comm); 

			//Now divide into columns
			int col_root=0,col_rank;
			//Divide the total list of processes into rows
			MPI_Comm  col_comm;
			MPI_Comm_split(outer_comm, key,color, &col_comm);
			MPI_Comm_rank(col_comm,&col_rank);

			if( row_rank== row_root  ) //ATTENTION IT IS ASSUMED THAT recv is initialize
				MPI_Reduce(&trecv[0],&recv[0],RecvSize, MPI_DOUBLE, MPI_SUM ,col_root,col_comm);

			MPI_Comm_free(&col_comm);
			MPI_Comm_free(&row_comm);
	
			if ( trecv != NULL ) delete[] trecv; trecv=NULL;
		}


int GetTasksPerProc(
	const int TotalNumTaks,
	MPI_Comm& comm )
{
	int ProcsInComm;
	MPI_Comm_size(comm, &ProcsInComm);

	const int 
	TaskPerProc = (TotalNumTaks+ProcsInComm-1)/ProcsInComm ;
	return TaskPerProc;
};


bool IsParallelizable(const int NumProc, const int ReqTasks, const double tol)
{	//tol value between 0 and 1, 0= perfect parallelizable, 1 worst case
	assert( NumProc>= ReqTasks);
	if ( NumProc%ReqTasks ==0 ) return true;
 
	return false;
};

size_t AssignRowColor(const int rank,const int size, const int ColorNum)
{
	const size_t ColorSize = (size+ColorNum-1)/ColorNum;
	return (rank/ColorSize+ColorSize)%ColorSize;
};

size_t AssignColColor(const int rank,const int size, const int ColorNum)
{
	return rank%ColorNum;
};

void GetParallelCommunicators(	MPI_Comm& oper_comm,
								MPI_Comm& rand_comm,
								MPI_Comm& ener_comm,
								const string configFilename,
								int world_rank,int  world_size
							 )
{

	if( world_rank==0)
		std::cout<<std::endl
				 <<"--------------PARALLELIZATION SCHEME--------------"<<std::endl
				 <<"The parallelzation scheme will be done in the "<<std::endl
				 <<"following way:"<<std::endl<<std::endl
				 <<"(1) Check if the parallelization by energy is enough "<<std::endl
				 <<"by checking if it consume of all the processes. "<<std::endl
				 <<"If it is, parallelize in energy and no further parallelization is performed."<<std::endl
				 <<"(1.1) If the parallelization in energy does not consume all processes "<<std::endl
				 <<"check if the parallelization by random vector is enough. If it is, then paralelize in random vector  parallelization is performed. "<<std::endl
				 <<"(1.2) If the number of cores is not an exact multiple of "<<std::endl
				 <<"the number of operators, then no parallelization is performed. "<<std::endl
				 <<"(1.3) If the parallelization is performed, the colors are defined"<<std::endl
				 <<"by the operators"<<std::endl
				 <<"(2) The second parallelization will be by random vector."<<std::endl
				 <<"which willl be done efficiently by making the total number of vectors variable."<<std::endl
				 <<"(2.1) If the number of random vector is larger that"
				 <<"the number of processor no further parallelization is performed."<<std::endl
				 <<"therefore the number of random vector can change."<<std::endl;

/*	
	std::ifstream configFile( configFilename.c_str() );

	//Get the necessary parameters to determine the most efficient 
	//parallelization
	std::string sim_label="";
	GetValue(configFile,"Label",sim_label);

	int NumRanVec=1;
	GetValue(configFile,"SamplingNumber", NumRanVec);

	double xmin,xmax,dx;
	GetValue(configFile,"Emin", xmin); 
	GetValue(configFile,"Emax", xmax);  
	GetValue(configFile,"dE", dx);  
	int NumEnergy=(xmax-xmin+dx)/dx;
*/
				 
	std::string sim_label="";
	size_t NumOper=1; 
	size_t NumRand=4;
	size_t NumEner=1;
	size_t TotNumTask=  NumOper*NumRand*NumEner;
	
	//DEFAULT BEHAVIOR: ALL COMMUNICATORS ARE EQUAL TO WORLD
	MPI_Comm_split(MPI_COMM_WORLD, world_rank, 0, &oper_comm);
	MPI_Comm_split(MPI_COMM_WORLD, world_rank, 0, &rand_comm);
	MPI_Comm_split(MPI_COMM_WORLD, world_rank, 0, &ener_comm);
	
	if( TotNumTask < world_size )
	{
		if( world_rank== 0)
			std::cout<<"\n--------------------WARNING--------------------"<<std::endl
					 <<"The total number of task is: "<<TotNumTask<<" "<<std::endl
					 <<"while the requested number of processes is: "<<world_size<<" "<<std::endl
					 <<"therefore no parallelization was performed"<<std::endl
					 <<"ALL PROCESSES WILL REPEAT THE SAME TASK"<<std::endl;
	return ;
	}


//----------------------PARALLLELIZATION OVER OPERATORS----------------//
	int ProcPerOp=-1; // -1 means SERIAL_CONF	
	if( NumOper > world_size)
	{
		if(world_rank==0) 
			std::cout<<"\n--------------------WARNING--------------------"<<std::endl
					 <<"Number of operators larger than"
					 <<"the number of processes."<<std::endl
					 <<"NO PARALLELIZATION PERFORMED "<<std::endl
					 <<"ALL PROCESSES WILL PERFORM THE SAME TASK"<<std::endl;
		MPI_Comm_split(MPI_COMM_WORLD, world_rank, 0, &oper_comm);
	}
	else
	{
		int rank=world_rank, size=world_size, tasksNum=NumOper;
		if( IsParallelizable(size,tasksNum,0) )
		{
			ProcPerOp=size/tasksNum;
			if(world_rank==0)
				std::cout	<<"\n       PARALLELIATION OVER VECTOR     "<<std::endl
							<<"Setting the parallelization of operators."<<std::endl
							<<"Using "<<ProcPerOp
							<<" process per group"<<std::endl;		
			//It will be assing to each operator a given set of processes
			int oper_color= AssignRowColor(rank,size,tasksNum); 
			MPI_Comm_split(MPI_COMM_WORLD, oper_color, rank, &oper_comm);
		}
	}
	if( ProcPerOp ==1)	//When all slots consumed, no further parallelization
	{	
		if( world_rank==0)			
		std::cout<<"\n The random vector parallelization consumed all slots."<<std::endl
				 <<"no further parallelization performes"<<std::endl;		
 		return ;
	}
	
	int oper_rank, oper_size;
	MPI_Comm_rank(oper_comm, &oper_rank);
	MPI_Comm_size(oper_comm, &oper_size);

	int oidx=oper_rank;
	if( oper_rank==1)
	do
	{
		std::cout<<"oidx: "<<oidx<<std::endl;
		oidx++;
	}
	while ( oper_size==1 && oidx< NumOper );
/*	
	//Parallelizing by operator
	
	int ProcPerOp=1;
	double tol=0; //The division of processes among operators has to be exact
	
	if( IsParallelizable(size,NumCol,tol) )
	{
		ProcPerOp=size/NumCol;
		if(world_rank==0)
			std::cout	<<"Setting the parallelization of operators."<<std::endl
						<<"Using "<<ProcPerOp
						<<" process per group"<<std::endl;		
		//It will be assing to each operator a given set of processes
		int oper_color= AssignRowColor(rank,size,NumCol); 
		MPI_Comm_split(MPI_COMM_WORLD, oper_color, rank, &oper_comm);

		if( ProcPerOp ==1)	//When all slots consumed, no further parallelization
		{}

	}
	else 
		if(world_rank==0) std::cout<<"\nNumber of processor not multiple"<<std::endl
								   <<"of the number of processes."<<std::endl
								   <<"NO PARALLELIZATION PERFORMED "<<std::endl
								   <<"COMM_RAND== COMM_WORLD"<<std::endl;
	
	int oper_rank, oper_size;
	MPI_Comm_rank(oper_comm, &oper_rank);
	MPI_Comm_size(oper_comm, &oper_size);
//	printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
//		world_rank, world_size, oper_rank, oper_size);


	
	if( ProcPerOp < NumColRand ) //No parallelization of random vector  is performed
	{
	}
	
	
	//Parallelizing by random_vector
	rank=oper_rank, size=oper_size, NumCol=NumColRand;
 	NumColRand = ( (NumColRand+ProcPerOp-1)/ProcPerOp)*ProcPerOp;
	int ProcPerRandPerOp=1;
	if( IsParallelizable(ProcPerOp,NumCol,0) )
	{
		ProcPerRandPerOp= NumCol/ProcPerOp;
		if(world_rank==0)
			std::cout	<<"Setting the parallelization of random vector."<<std::endl
						<<"Using "<<ProcPerRandPerOp
						<<" process per group"<<std::endl;		
		MPI_Comm_dup(MPI_COMM_WORLD,&rand_comm);
	}
	else 
		if(world_rank==0) std::cout<<"\nNumber of processor not multiple"<<std::endl
								   <<"of the number of random vector ."<<std::endl
								   <<"NO PARALLELIZATION PERFORMED "<<std::endl
								   <<"COMM_RAND== COMM_WORLD"<<std::endl;

	int rand_rank, rand_size;
	MPI_Comm_rank(rand_comm, &rand_rank);
	MPI_Comm_size(rand_comm, &rand_size);0

	if( ProcPerRandPerOp < NumColEner ) //No parallelization of energy is performed
	{
	}
	
	
	//Parallelizing by energy
	rank=ener_rank, size=ener_size, NumCol=NumColEner;
	int ProcEnerPerRandPerOp=1;

	ProcEnerPerRandPerOp
		ProcPerRandPerOp= NumCol/ProcPerOp;
		if(world_rank==0)
			std::cout	<<"Setting the parallelization of random vector."<<std::endl
						<<"Using "<<ProcPerRandPerOp
						<<" process per group"<<std::endl;		
		MPI_Comm_dup(MPI_COMM_WORLD,&rand_comm);
	}
	else 
		if(world_rank==0) std::cout<<"\nNumber of processor not multiple"<<std::endl
								   <<"of the number of random vector ."<<std::endl
								   <<"NO PARALLELIZATION PERFORMED "<<std::endl
								   <<"COMM_RAND== COMM_WORLD"<<std::endl;

	int rand_rank, rand_size;
	MPI_Comm_rank(rand_comm, &rand_rank);
	MPI_Comm_size(rand_comm, &rand_size);0
	
*/

/*
	int oper_rank,oper_size;
	MPI_Comm_rank(oper_comm, &oper_rank);
	MPI_Comm_size(oper_comm, &oper_size);
	printf("WORLD RANK/SIZE: %d/%d \t MY OPER_RANK IS=%d OUT OF %d\n",
		world_rank, world_size, oper_rank, oper_size);	

//	if(oper_rank==0)
//		std::cout<<"I am rank "<<oper_rank<<"/"<<oper_size<<" of color "<<oper_color<<" and global rank "<<world_rank<<std::endl;

	const int OperPerCol=(NumColOper+oper_size-1)/oper_size;

	for( int op=oper_rank;op<(oper_rank+1)*OperPerCol; op++)
//	if( oper_rank==0)
	{
	//	std::cout<<"The OpColor group "<<oper_rank<<"  will compute "<<OperPerCol<<" out of the total  "<<NumColOper<<" operators"<<std::endl;
	}*/
/*	




	if( oper_rank== 0)
		std::cout<<"I have color= "<<oper_rank<<std::endl;*/

	// Get the rank and size in the original communicator
//	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


//	int color = world_rank / 2; // Determine color based on row

	// Split the communicator based on the color and use the
	// original rank for ordering
//	MPI_Comm row_comm;
//	MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

//	int row_rank, row_size;
//	MPI_Comm_rank(row_comm, &row_rank);
//	MPI_Comm_size(row_comm, &row_size);

//	printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
//		world_rank, world_size, row_rank, row_size);

//	MPI_Comm_free(&row_comm);

/*

	//Get the number of task for each type of parallelization
	const int 
	world_root=0,
	rand_taskPerProc = ( world_size +NumRanVec -1 )/ world_size ,
	ener_taskPerProc = ( world_size +NumEnergy -1 )/ world_size ;
	
	
	if( world_rank== world_root)
	{
		std::cout<<"The parallelization by random vector"<<std::endl
				 <<" leads to a number of task/proc of "<<rand_taskPerProc<<std::endl
				 <<"The parallelization by energy"<<std::endl
				 <<" leads to a number of task/proc of "<<ener_taskPerProc<<std::endl;
		if ( ener_taskPerProc>=rand_taskPerProc)
			std::cout<<"Using energy parallelization (Preferred)"<<std::endl;
		else
		{
			int newNumRanVec=rand_taskPerProc*world_size;
			std::cout<<"Using random vector parallelization"<<std::endl;
			std::cout<<"The number of random vectors for the simulation "
				 <<sim_label<<std::endl
				 <<" was changed from "<<NumRanVec<<" to "<<newNumRanVec
				 <<std::endl<<" in order to maximize the efficiency"
				 <<" of core division"<<std::endl;
			NumRanVec=newNumRanVec;
		}
	}

	int rand_color=world_rank, ener_color=world_rank;
	if ( ener_taskPerProc >= rand_taskPerProc) ener_color=1;
	else  rand_color =1;

	MPI_Comm_split(MPI_COMM_WORLD, rand_color, world_rank, &rand_comm);
	MPI_Comm_split(MPI_COMM_WORLD, ener_color, world_rank, &ener_comm);
*/
};

#endif
