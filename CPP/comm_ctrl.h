#ifndef COMM_CTRL
#define COMM_CTRL

#include <stdlib.h>
#include <mpi.h>
#include "decomposition.h"

class comm_ctrl{
public: 
	comm_ctrl(decomposition *d, double *rhs, double *u, double *U_up, 
						int bottomrank, int toprank, int leftrank, int rightrank,
						int mygrp, int bottomgrp, int topgrp, int leftgrp, int rightgrp);
	~comm_ctrl();
	void send();
	void receive();
	void wait();
	void compile_solution();
private:
	int myRank, bottomRank, topRank, leftRank, rightRank,
			myGrp, bottomGrp, topGrp, leftGrp, rightGrp;
	decomposition *D;
	double *RHS, *U, *u_up;
//
	void compile_updates();
	void compile_update_toBottom();
	void compile_update_toTop();
	void compile_update_toLeft();
	void compile_update_toRight();
	void send_updates();
	void send_update_toBottom();
	void send_update_toTop();
	void send_update_toLeft();
	void send_update_toRight();
	void receive_updates();
	void receive_update_toBottom();
	void receive_update_toTop();
	void receive_update_toLeft();
	void receive_update_toRight();
};

#endif