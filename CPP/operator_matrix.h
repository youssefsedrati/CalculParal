#ifndef OPERATOR_MATRIX_H
#define OPERATOR_MATRIX_H

class operator_matrix{
public:
	operator_matrix(int Nx, int Ny, double Lx, double Ly, double D);
	operator_matrix();
	~operator_matrix();
	double Aii();
	double Cx();
	double Cy();
private:
	double _Aii, _Cx, _Cy, dx, dy;
};

#endif