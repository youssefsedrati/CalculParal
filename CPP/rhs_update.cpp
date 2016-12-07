#include "operator_matrix.h"
#include "decomposition.h"
#include "rhs_update.h"
using namespace std;

RHS_Update::RHS_update(decomposition *d,operator_matrix *a,double *rhs) {	
	D = d; A = a; RHS = rhs;
}

RHS_update::~RHS_update(){
}


void RHS_Update::RHS_update_top(double *u_up){
	index = D->get_index_global_top();
	for(int i=0;i<D->get_myNx();i++){ 
      RHS[index[i]] =RHS[index[i]]-u_up[i]*A->Cy();
  }
}

void RHS_Update::RHS_update_bottom(double *u_up){
	index = D->get_index_global_bottom();
	for(int i=0;i<D->get_myNx();i++){ 
      RHS[index[i]] =RHS[index[i]]-u_up[i]*A->Cy();
  }
}

void RHS_Update::RHS_update_right(double *u_up){
	index = D->get_index_global_right();
	for(int i=0;i<D->get_myNy();i++){ 
      RHS[index[i]] =RHS[index[i]]-u_up[i]*A->Cx();
  }
}

void RHS_Update::RHS_update_left(double *u_up){
	index = D->get_index_global_left();
	for(int i=0;i<D->get_myNy();i++){ 
      RHS[index[i]] =RHS[index[i]]-u_up[i]*A->Cx();
  }
}

