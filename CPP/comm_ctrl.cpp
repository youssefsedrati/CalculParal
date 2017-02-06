#include "comm_ctrl.h"
#include "decomposition.h"
#include "operator_matrix.h"

using namespace std; 

/* comm_ctrl holds information of adjacency for each process
   in directions (top,bottom,left,right);
	 according to which updates are sent (boundary values)
	 and received (incorporate these values into RHS);
	 
	 u_top/bot... are vectors used in transmissions;
*/
comm_ctrl::comm_ctrl(decomposition *d, operator_matrix *a,
			double *rhs, double *rhs_up){
	D=d; A=a;RHS=rhs; RHS_up=rhs_up;
	myN = D->get_myN(); myNx = D->get_myNx(); myNy = D->get_myNy(); 
	u_bot = (double*)malloc(myNx*sizeof(double));
	u_top = (double*)malloc(myNx*sizeof(double));
	u_left = (double*)malloc(myNy*sizeof(double));
	u_right = (double*)malloc(myNy*sizeof(double));
	init_neighbour_ranks();
	init_group_behaviour();
}

comm_ctrl::~comm_ctrl(){
	cleanup();
}

void comm_ctrl::cleanup(){
	if(u_bot) free(u_bot);
	if(u_top) free(u_top);
	if(u_left) free(u_left);
	if(u_right) free(u_right);
}

void comm_ctrl::send(double *U){
	send_updates(U);
}

void comm_ctrl::receive(){
	if(is_myGroup_waiting)
		receive_updates();
	else
		is_myGroup_waiting=true;
}

// fills solution vector at end of calculus
void comm_ctrl::compile_solution(double *U){
	MPI_Allreduce(MPI_IN_PLACE, U, D->get_N(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

// sum of all squared norms
void comm_ctrl::cumulate_dist_squared(double *dsq){
	MPI_Allreduce(MPI_IN_PLACE, dsq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void comm_ctrl::init_neighbour_ranks(){
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	bottomRank= (D->get_myRank_y()>0) 									? myRank - D->get_N_procs_x() : -1;
	topRank= (D->get_myRank_y()<D->get_N_procs_y()-1) 	? myRank + D->get_N_procs_x() : -1;
	leftRank= (D->get_myRank_x()>0) 										? myRank-1 : -1;
	rightRank= (D->get_myRank_x()<D->get_N_procs_x()-1) ? myRank+1 : -1;
}

void comm_ctrl::init_group_behaviour(){
	am_I_in_rootGroup = ( (D->get_myRank_x()+D->get_myRank_y())%2==0 ) ? true:false;
	is_myGroup_waiting = !am_I_in_rootGroup;
}

// sends integer to sendRank who receives from his recvRank (this process)
void comm_ctrl::comm_neighbour_int(int sendRank,int recvRank,int* msg,int* recv_buf, int label){
	if(sendRank>=0){
		MPI_Send(msg,1,MPI_INTEGER,sendRank,label,MPI_COMM_WORLD);
		//cout << myRank << ". " << *msg << " sent to " << sendRank << " as " << recvRank <<endl;
	}
	if(recvRank>=0){
		MPI_Recv(recv_buf,1,MPI_INTEGER,recvRank,label,MPI_COMM_WORLD,&mpi_stat);
		//cout << myRank << ". " << *recv_buf << " recvd from " << recvRank <<endl;
	}
}

/* main send routine
	 uses MPI_Send to transmit Dirichlet boundary data
*/
void comm_ctrl::send_updates(double *U){
	send_update_toBottom(U);
	send_update_toTop(U);
	send_update_toLeft(U);
	send_update_toRight(U);
}

void comm_ctrl::send_update_toBottom(double *U){
	if(bottomRank<0) return;
	int *idx = D->get_index_global_bottom();
	for(int i=0;i<myNx;++i){
		int j = idx[i];
		u_bot[i] = U[j];
	}
	MPI_Send(u_bot,myNx,MPI_DOUBLE,bottomRank,100,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toTop(double *U){
	if(topRank<0) return;
	int *idx = D->get_index_global_top();
	for(int i=0;i<myNx;++i){
		int j = idx[i];
		u_top[i] = U[j];
	}
	MPI_Send(u_top,myNx,MPI_DOUBLE,topRank,200,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toLeft(double *U){
	if(leftRank<0) return;
	int *idx = D->get_index_global_left();
	for(int i=0;i<myNy;++i){
		int j = idx[i];
		u_left[i] = U[j];
	}
	MPI_Send(u_left,myNy,MPI_DOUBLE,leftRank,300,MPI_COMM_WORLD);
}

void comm_ctrl::send_update_toRight(double *U){
	if(rightRank<0) return;
	int *idx = D->get_index_global_right();
	for(int i=0;i<myNy;++i){
		int j = idx[i];
		u_right[i] = U[j];
	}
	MPI_Send(u_right,myNy,MPI_DOUBLE,rightRank,400,MPI_COMM_WORLD);
}

/* main receive routine 
	 uses MPI_Recv to obtain Dirichlet boundary data
	 data is then worked into the RHS update vector.
*/
void comm_ctrl::receive_updates(){
	receive_update_fromBottom();
	receive_update_fromTop();
	receive_update_fromLeft();
	receive_update_fromRight();
}

void comm_ctrl::receive_update_fromBottom(){
	if(bottomRank<0) return;
	MPI_Recv(u_bot, myNx, MPI_DOUBLE, bottomRank, 200, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_bottom();
	for(int i=0;i<myNx;++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cy() * u_bot[i];
	}
}

void comm_ctrl::receive_update_fromTop(){
	if(topRank<0) return;
	MPI_Recv(u_top,myNx, MPI_DOUBLE, topRank, 100, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_top();
	for(int i=0;i<myNx;++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cy() * u_top[i];
	}
}

void comm_ctrl::receive_update_fromLeft(){
	if(leftRank<0) return;
	MPI_Recv(u_left,myNy, MPI_DOUBLE, leftRank, 400, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_left();
	for(int i=0;i<myNy;++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cx() * u_left[i];
	}
}

void comm_ctrl::receive_update_fromRight(){
	if(rightRank<0) return;
	MPI_Recv(u_right,myNy, MPI_DOUBLE, rightRank, 300, MPI_COMM_WORLD, &mpi_stat);
	int *idx = D->get_index_global_right();
	for(int i=0;i<myNy;++i){
		int j = idx[i];
		RHS_up[j] = RHS[j] - A->Cx() * u_right[i];
	}
}