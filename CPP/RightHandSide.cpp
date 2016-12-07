#include <iostream>
#include <fstream>
#include <string>
#include "operator_matrix.h"
#include "decomposition.h"
#include "RightHandSide.h"
#include "tools.h"
#include "cdt_bords.h"
using namespace std;

RightHandSide::RightHandSide() {	

	double (*f)(double, double, double) = cdt[cdt_choisie].f;
  	double (*g)(double, double, double) = cdt[cdt_choisie].g;
  	double (*h)(double, double, double) = cdt[cdt_choisie].h;

}

RightHandSide::~RightHandSide(){

}


void RightHandSide::initRHS(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS)
{

int i,k,j,posx,posy;
 for( i = 1; i<=N; i++ ){
    nloc(&j,&k,i,Nx);
    posx = k*dx;
    posy = j*dy;
    RHS[i] = f(posx,posy,0.0);
  }

}


void RightHandSide::computeRHS_alternate_update_top(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS){

int i,k,j,posx,posy;
	for(i = N-Nx+1;i<=N;i++ ){ 
      nloc(&j,&k,i,Nx);
      posx = k*dx;
      posy = j*dy;
      RHS[i] =RHS[i]-g(posx,1.0,0.0)*Cy;
    }

}


void RightHandSide::computeRHS_alternate_update_bottom(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS){

int i,k,j,posx,posy;	
	 for( i = 1; i<= Nx; i++ ){
  		nloc(&j,&k,i,Nx);
  		posx = k*dx;
  		posy = j*dy;
  		RHS[i] = RHS[i]-g(posx,0.0,0.0)*Cy;
      }

}
void RightHandSide::computeRHS_alternate_update_right(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS){
	// Ligne du bas
	RHS[1]  = RHS[1]  -h(0.0,dy,0.0)*Cx;
	// Ligne du haut
	RHS[N-Nx+1] = RHS[N-Nx+1] -h(0.0,k*dy,0.0)*Cx;
}
void RightHandSide::computeRHS_alternate_update_left(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS){
	// Ligne du bas.
	RHS[Nx] = = RHS[Nx] -h(1.0,dy,0.0)*Cx;
	// Ligne du haut
	RHS[N] = RHS[N] -h(1.0,k*dy,0.0)*Cx;

}

void computeRHS_alternate_update_middle(int N, int Nx, double dx, double dy, double Cx, double Cy, double *RHS)
{

	int j = 1 + Nx;
	int M = N/Nx;
	for( i = 2; i<= M-1; i++ ){
    	nloc(&k,&l,j,Nx);
    	RHS[j] = RHS[j] -h(0.0,k*dy,0.0)*Cx;
    	RHS[j+Nx-1] = RHS[j+Nx-1] -h(1.0,k*dy,0.0)*Cx;
    	j = 1 + (i)*Nx;
  	}
}
