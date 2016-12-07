#ifndef RHS_H
#define RHS_H

#include <iostream>
#include <fstream>
#include <string>
#include "operator_matrix.h"
#include "decomposition.h"
#include "tools.h"
#include "cdt_bords.h"

class RightHandSide{

public:
	RightHandSide();

	~RightHandSide();

private:
	int cdt_choisie= 1;
	int N , Nx,M;
	double dx,dy, Cx,Cy;
	double *RHS;

	void initRHS(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS);
	void computeRHS_alternate_update_top(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS);
	void computeRHS_alternate_update_bottom(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS);
	void computeRHS_alternate_update_right(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS);
	void computeRHS_alternate_update_left(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS);
	void computeRHS_alternate_update_middle(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS);



};








#endif