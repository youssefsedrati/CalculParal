#include "comm_ctrl.h"
#include "decomposition.h"

comm_ctrl::comm_ctrl(decomposition *d, double *rhs, double *u, double *U_up, 
			int bottomrank, int toprank, int leftrank, int rightrank,
			int mygrp, int bottomgrp, int topgrp, int leftgrp, int rightgrp){				
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	D = d; u_up = U_up; RHS = rhs, U = u;
	bottomRank=bottomrank; topRank=toprank; leftRank=leftrank; rightRank=rightrank;
	myGrp=mygrp; bottomGrp=bottomgrp; topGrp=topgrp; leftGrp=leftgrp; rightGrp=rightgrp;
}

comm_ctrl::~comm_ctrl(){}

void comm_ctrl::send(){
	
}

void comm_ctrl::receive(){
	
}

void comm_ctrl::wait(){
	
}

void comm_ctrl::compile_solution(){
	
}

void comm_ctrl::compile_updates(){
	
}

void comm_ctrl::compile_update_toBottom(){
	
}

void comm_ctrl::compile_update_toTop(){
	
}

void comm_ctrl::compile_update_toLeft(){
	
}

void comm_ctrl::compile_update_toRight(){
	
}

void comm_ctrl::send_updates(){
	
}

void comm_ctrl::send_update_toBottom(){
	
}

void comm_ctrl::send_update_toTop(){
	
}

void comm_ctrl::send_update_toLeft(){
	
}

void comm_ctrl::send_update_toRight(){
	
}

void comm_ctrl::receive_updates(){
	
}

void comm_ctrl::receive_update_Bottom(){
	
}

void comm_ctrl::receive_update_Top(){
	
}

void comm_ctrl::receive_update_Left(){
	
}

void comm_ctrl::receive_update_Right(){
	
}