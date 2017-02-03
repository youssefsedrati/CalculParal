#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>
#include <string>
#include "operator_matrix.h"
#include "decomposition.h"

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define UNUSED(x) (void)(x)

/** Calcule le i et j correspondant au n passé en paramètre, et les affecte
*   aux pointeurs i et j.
*/
void nloc(int *i, int *j, int n, int Nx);



/** Calcule le second membre en fonction des conditions
  * aux bords définies globalement 
  */
void RightHandSide(int N, int Nx, int M, double dx, double dy,
				   double Cx, double Cy,double *RHS);


/* Produit matrice vecteur pour une matrice tridiagonale par bloc
! Produit matrice vecteur dans le cas ou A pentadiagonale de la forme:
!
! A = B C             matrice pentadiagonale (m*Nx,m*Nx)
!     C B C
!       C B C
!         . . .
!           . . .
!             C B
! avec 
! C = Cy*Id            matrice diagonale (Nx,Nx)
! 
! B = Aii Cx           matrice tridiagonale (Nx,Nx)
!     Cx  Aii Cx
!          .   .   .
!              Cx Aii
*/
//void matvec(double Aii,double Cx,double Cy,int Nx,int Ny,double *Uold,double *U);
void matvec(operator_matrix a , decomposition *D, double *Uold,double *U);
/** returns the square l2-distance of two parsed vectors of length N+1
	*/
double dist_squared(int N,double* U,double* V);

#endif
