#include "operator_matrix.h"

operator_matrix::operator_matrix(int Nx, int Ny, double Lx, double Ly, double D){
	dx  = Lx/(1.0 + Nx);
	dy  = Ly/(1.0 + Ny);
	_Aii = 2.0*D/(dx*dx)+ 2.0*D/(dy*dy);
	_Cx  = -1.0*D/(dx*dx);
	_Cy  = -1.0*D/(dy*dy);
}

operator_matrix::operator_matrix(){
	
}

operator_matrix::~operator_matrix(){
	
}

double operator_matrix::Aii(){
	return _Aii;
}

double operator_matrix::Cx(){
	return _Cx;
}

double operator_matrix::Cy(){
	return _Cy;
}