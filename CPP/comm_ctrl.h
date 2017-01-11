#ifndef COMM_CTRL
#define COMM_CTRL

#include <stdlib.h>
#include <mpi.h>
#include "decomposition.h"

class comm_ctrl{
public: 
	comm_ctrl(decomposition *d, double *rhs, double *rhs_up, double *u, double *U_up, 
						int bottomrank, int toprank, int leftrank, int rightrank);
	~comm_ctrl();
	void send();
	void receive();
	void compile_solution();
private:
	int myRank, bottomRank, topRank, leftRank, rightRank;
	decomposition *D;
	double *RHS, *RHS_up, *U, *u_up;
	MPI_Status mpi_stat;
//
	void send_updates();
	void send_update_toBottom();
	void send_update_toTop();
	void send_update_toLeft();
	void send_update_toRight();
	void receive_updates();
	void receive_update_fromBottom();
	void receive_update_fromTop();
	void receive_update_fromLeft();
	void receive_update_fromRight();
};

#endif