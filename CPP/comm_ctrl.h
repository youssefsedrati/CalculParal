#ifndef COMM_CTRL
#define COMM_CTRL


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <mpi.h>

#include "decomposition.h"
#include "operator_matrix.h"

class comm_ctrl{
public: 
	comm_ctrl(decomposition *d, operator_matrix *a,
						double *rhs, double *rhs_up, double *U_up);
	~comm_ctrl();
	void send(double *U);
	void receive();
	void compile_solution(double *U);
	void cumulate_dist_squared(double* dsq);
private:
	int myRank, bottomRank, topRank, leftRank, rightRank;
	bool am_I_in_rootGroup, is_myGroup_waiting;
	decomposition *D;
	operator_matrix *A;
	double *RHS, *RHS_up, *u_up;
	MPI_Status mpi_stat;
//
	void init_neighbour_ranks();
	void init_group_behaviour();
	void send_updates(double *u);
	void send_update_toBottom(double *u);
	void send_update_toTop(double *u);
	void send_update_toLeft(double *u);
	void send_update_toRight(double *u);
	void receive_updates();
	void receive_update_fromBottom();
	void receive_update_fromTop();
	void receive_update_fromLeft();
	void receive_update_fromRight();
};

#endif