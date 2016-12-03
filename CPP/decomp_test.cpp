#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "decomposition.h"

using namespace std;

int main(){
	int myrank=4, nb_procs=8, nb_procs_x=4, nx=10, ny=10;
	for(int j=0;j<nb_procs;++j){
		decomposition D(j, nb_procs, nb_procs_x, nx, ny);
		int* index = D.get_index_global_right();
		cout << "#" << j << ". " << D.get_myNx()<< ": ";
		for(int i=0;i<D.get_myNy();++i)
			cout << index[i] << " ";
		cout << endl;
	}
	return 0;
}