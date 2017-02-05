#ifndef OPERATOR_MATRIX_H
#define OPERATOR_MATRIX_H

class operator_matrix{
public:
	operator_matrix(int Nx, int Ny, double Lx, double Ly, double D, bool NeumannBC);
	operator_matrix();
	~operator_matrix();
	double Aii();
	double Cx();
	double Cy();
	double dx();
	double dy();
private:
	double _Aii, _Cx, _Cy, _dx, _dy;
};

#endif