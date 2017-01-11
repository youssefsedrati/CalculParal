#include "comm_ctrl.h"
#include "decomposition.h"

comm_ctrl::comm_ctrl(decomposition *d, double *rhs, double *rhs_up, double *u, double *U_up, 
			int bottomrank, int toprank, int leftrank, int rightrank){				
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	D = d; u_up = U_up; RHS = rhs, U = u; RHS_up = rhs_up;
	bottomRank=bottomrank; topRank=toprank; leftRank=leftrank; rightRank=rightrank;
}

comm_ctrl::~comm_ctrl(){}

void comm_ctrl::send(){
	send_updates();
}

void comm_ctrl::receive(){
	receive_updates();
}

void comm_ctrl::compile_solution(){
	MPI_Allreduce(MPI_IN_PLACE, U, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void comm_ctrl::send_updates(){
	send_update_toBottom();
	send_update_toTop();
	send_update_toLeft();
	send_update_toRight();
}

void comm_ctrl::send_update_toBottom(){
	if(bottomRank<0) return;
	int *index = D->get_index_global_top();
	for(size_t i=0;i<D->get_myNx();++i)
		u_up[index[i]] = U[index[i]];
	MPI_Send(u_up,D->get_myN(),MPI_DOUBLE,bottomRank,100,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toTop(){
	if(topRank<0) return;
	int *index = D->get_index_global_bottom();
	for(size_t i=0;i<D->get_myNx();++i)
		u_up[index[i]] = U[index[i]];
	MPI_Send(u_up,D->get_myN(),MPI_DOUBLE,topRank,200,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toLeft(){
	if(leftRank<0) return;
	int *index = D->get_index_global_right();
	for(size_t i=0;i<D->get_myNy();++i)
		u_up[index[i]] = U[index[i]];
	MPI_Send(u_up,D->get_myN(),MPI_DOUBLE,leftRank,300,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toRight(){
	if(rightRank<0) return;
	int *index = D->get_index_global_left();
	for(size_t i=0;i<D->get_myNy();++i)
		u_up[index[i]] = U[index[i]];
	MPI_Send(u_up,D->get_myN(),MPI_DOUBLE,rightRank,400,MPI_COMM_WORLD);
}

void comm_ctrl::receive_updates(){
	receive_update_fromBottom();
	receive_update_fromTop();
	receive_update_fromLeft();
	receive_update_fromRight();
}

void comm_ctrl::receive_update_fromBottom(){
	if(bottomRank<0) return;
	MPI_Recv(u_up, D->get_myN(), MPI_DOUBLE, bottomRank, 100, MPI_COMM_WORLD, &mpi_stat);
	int *index = D->get_index_global_bottom();
	for(size_t i=0;i<D->get_myNx();++i){
		RHS_up[index[i]] = RHS[index[i]]
			+D->get_myN()*D->get_myN() * u_up[index[i]-D->get_myNx()]; 
	}
}

void comm_ctrl::receive_update_fromTop(){
	if(topRank<0) return;
	MPI_Recv(u_up, D->get_myN(), MPI_DOUBLE, topRank, 200, MPI_COMM_WORLD, &mpi_stat);
	int *index = D->get_index_global_top();
	for(size_t i=0;i<D->get_myNx();++i){
		RHS_up[index[i]] = RHS[index[i]]
			+D->get_myN()*D->get_myN() * u_up[index[i]+D->get_myNx()]; 
	}
}

void comm_ctrl::receive_update_fromLeft(){
	if(leftRank<0) return;
	MPI_Recv(u_up, D->get_myN(), MPI_DOUBLE, leftRank, 300, MPI_COMM_WORLD, &mpi_stat);
	int *index = D->get_index_global_left();
	for(size_t i=0;i<D->get_myNy();++i){
		RHS_up[index[i]] = RHS[index[i]]
			+D->get_myN()*D->get_myN() * u_up[index[i]-1]; 
	}
}

void comm_ctrl::receive_update_fromRight(){
	if(rightRank<0) return;
	MPI_Recv(u_up, D->get_myN(), MPI_DOUBLE, rightRank, 400, MPI_COMM_WORLD, &mpi_stat);
	int *index = D->get_index_global_bottom();
	for(size_t i=0;i<D->get_myNy();++i){
		RHS_up[index[i]] = RHS[index[i]]
			+D->get_myN()*D->get_myN() * u_up[index[i]+1]; 
	}
}