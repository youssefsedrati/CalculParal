#include "operator_matrix.h"

/* operator_matrix holds the contents of the discrete Laplacian;
	 three values in total: diagonal, close side diagonals, far side diagonals
*/
operator_matrix::operator_matrix(int Nx, int Ny, double Lx, double Ly, double D, bool NeumannBC){
	if(NeumannBC){
		_dx  = Lx/(-1 + Nx);
		_dy  = Ly/(-1 + Ny);
	}else{
		_dx  = Lx/(1 + Nx);
		_dy  = Ly/(1 + Ny);
	}	
	_Aii = 2.0*D/(_dx*_dx)+ 2.0*D/(_dy*_dy);
	_Cx  = -1.0*D/(_dx*_dx);
	_Cy  = -1.0*D/(_dy*_dy);
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

double operator_matrix::dx(){
	return _dx;
}

double operator_matrix::dy(){
	return _dy;
}