#ifndef RHS_H_ 
#define RHS_H_ 
 
#include <math.h> 
#include "decomposition.h" 
#include "operator_matrix.h" 
 
using namespace std; 
 
double f(double x,double y){ 
  return 2*(x*(1-x) + y*(1-y) ); 
}  

double f1(double x,double y){ 
  return f(x,y)/4;
} 
 
double g(double x,double y){ 
  return sin(x)+cos(y); 
} 
 
double null(double x,double y){ 
  return 0*x*y; 
} 
 
double one(double x,double y){ 
  return 1+ 0*x*y; 
} 
 
void fill_RHS_force(decomposition *D,operator_matrix *A,  
      double *RHS,double(*func)(double,double)){ 
  int Nx=D->get_Nx(), Ny=D->get_Ny(); 
  for(int i=0;i<Ny;++i) 
    for(int j=0;j<Nx;++j){ 
      double x = (double)j/(Nx+1), y = (double)i/(Ny+1); 
      RHS[j+Nx*i] = func(x,y); 
    } 
} 
 
void fill_RHS_DirichletBC(decomposition *D,operator_matrix *A,  
      double *RHS,double(*func)(double,double)){ 
  int Nx=D->get_Nx(), Ny=D->get_Ny(); 
  double Cx=A->Cx(), Cy=A->Cy(); 
  for(int i=0;i<Nx;++i){ 
    double x = (double)i/(Nx-1), y = 0; 
    RHS[i] -= func(x,y)*Cx; 
    y = 1; 
    RHS[i+Nx*(Nx-1)]-= func(x,y)*Cx; 
  } 
  for(int i=0;i<Ny;++i){ 
    double x = 0, y =(double)i/(Ny-1); 
    RHS[i*Nx] -= func(x,y)*Cy; 
    x = 1; 
    RHS[(i+1)*Nx-1]-= func(x,y)*Cy; 
  } 
} 
 
void fill_RHS_NeumannBC(decomposition *D,operator_matrix *A,  
      double *RHS,double(*func)(double,double)){ 
  int Nx=D->get_Nx(), Ny=D->get_Ny(); 
  double dx=A->dx(), dy=A->dy(); 
  for(int i=1;i<Nx-1;++i){ 
    double x = (double)i/(Nx-1), y = 0; 
    RHS[i] -= func(x,y)/dx; 
    y = 1; 
    RHS[i+Nx*(Nx-1)]-= func(x,y)/dx; 
  } 
  for(int i=1;i<Ny-1;++i){ 
    double x = 0, y =(double)i/(Ny-1); 
    RHS[i*Nx] -= func(x,y)/dy; 
    x = 1; 
    RHS[(i+1)*Nx-1]-= func(x,y)/dy; 
  } 
} 
 
#endif