#include "comm_ctrl.h"
#include "decomposition.h"
#include "operator_matrix.h"

comm_ctrl::comm_ctrl(decomposition *d, operator_matrix *a,
			double *rhs, double *rhs_up, double *u, double *U_up){			
	D=d; A=a, u_up=U_up; RHS=rhs, U=u; RHS_up=rhs_up;
	init_neighbour_ranks();
	init_group_behaviour();
}

comm_ctrl::~comm_ctrl(){}

void comm_ctrl::send(){
	send_updates();
}

void comm_ctrl::receive(){
	if(is_myGroup_waiting)
		receive_updates();
	else
		is_myGroup_waiting=true;
}

void comm_ctrl::compile_solution(){
	MPI_Allreduce(MPI_IN_PLACE, U, D->get_N(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void comm_ctrl::init_neighbour_ranks(){	
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	bottomRank= (D->get_myRank_y()>0) 									? myRank - D->get_N_procs_x() : -1; 
	topRank= (D->get_myRank_y()<D->get_N_procs_y()-1) 	? myRank + D->get_N_procs_x() : -1; 
	leftRank= (D->get_myRank_x()>0) 										? myRank-1 : -1; 
	rightRank= (D->get_myRank_x()<D->get_N_procs_x()-1) ? myRank+1 : -1;
	//std::cout << "#" << myRank << ". v " << bottomRank << ". ^ " << topRank << ". < " << leftRank << ". > " << rightRank << ".\n";
}

void comm_ctrl::init_group_behaviour(){
	am_I_in_rootGroup = ( (D->get_myRank_x()+D->get_myRank_y())%2==0 ) ? true:false;
	is_myGroup_waiting = !am_I_in_rootGroup;
}

void comm_ctrl::send_updates(){
	send_update_toBottom();
	send_update_toTop();
	send_update_toLeft();
	send_update_toRight();
}

void comm_ctrl::send_update_toBottom(){
	if(bottomRank<0) return;
	int *idx = D->get_index_global_bottom();
	for(size_t i=0;i<D->get_myNx();++i){
		int j = idx[i];
		u_up[j] = U[j];
	}
	MPI_Send(u_up,D->get_N(),MPI_DOUBLE,bottomRank,100,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toTop(){
	if(topRank<0) return;
	int *idx = D->get_index_global_top();
	for(size_t i=0;i<D->get_myNx();++i){
		int j = idx[i];
		u_up[j] = U[j];
	}
	MPI_Send(u_up,D->get_N(),MPI_DOUBLE,topRank,200,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toLeft(){
	if(leftRank<0) return;
	int *idx = D->get_index_global_left();
	for(size_t i=0;i<D->get_myNy();++i){
		int j = idx[i];
		u_up[j] = U[j];
	}
	MPI_Send(u_up,D->get_N(),MPI_DOUBLE,leftRank,300,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toRight(){
	if(rightRank<0) return;
	int *idx = D->get_index_global_right();
	for(size_t i=0;i<D->get_myNy();++i){
		int j = idx[i];
		u_up[j] = U[j];
	}
	MPI_Send(u_up,D->get_N(),MPI_DOUBLE,rightRank,400,MPI_COMM_WORLD);
}

void comm_ctrl::receive_updates(){
	receive_update_fromBottom();
	receive_update_fromTop();
	receive_update_fromLeft();
	receive_update_fromRight();
}

void comm_ctrl::receive_update_fromBottom(){
	if(bottomRank<0) return;
	MPI_Recv(u_up, D->get_N(), MPI_DOUBLE, bottomRank, 200, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_bottom(),
			offset = D->get_Nx();
	for(size_t i=0;i<D->get_myNx();++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cy() * u_up[j - offset]; 
	}
}

void comm_ctrl::receive_update_fromTop(){
	if(topRank<0) return;
	MPI_Recv(u_up, D->get_N(), MPI_DOUBLE, topRank, 100, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_top(),
			offset = D->get_Nx();
	for(size_t i=0;i<D->get_myNx();++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cy() * u_up[j + offset]; 
	}
}

void comm_ctrl::receive_update_fromLeft(){
	if(leftRank<0) return;
	MPI_Recv(u_up, D->get_N(), MPI_DOUBLE, leftRank, 400, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_left();
	for(size_t i=0;i<D->get_myNy();++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cx() * u_up[j-1]; 
	}
}

void comm_ctrl::receive_update_fromRight(){
	if(rightRank<0) return;
	MPI_Recv(u_up, D->get_N(), MPI_DOUBLE, rightRank, 300, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_bottom();
	for(size_t i=0;i<D->get_myNy();++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cx() * u_up[j+1]; 
	}
}