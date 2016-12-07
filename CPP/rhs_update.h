#ifndef RHS_H
#define RHS_H

#include <iostream>
#include "operator_matrix.h"
#include "decomposition.h"

class RHS_update{

public:
	RHS_update(decomposition *d,operator_matrix *a,double *rhs);
	~RHS_update();

private:
	decomposition *D;
	operator_matrix *A;
	double *RHS;
	int *index;

	void RHS_update_top(double *u_up);
	void RHS_update_bottom(double *u_up);
	void RHS_update_right(double *u_up);
	void RHS_update_left(double *u_up);
};

#endif