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