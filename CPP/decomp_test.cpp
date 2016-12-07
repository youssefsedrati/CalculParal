#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "decomposition.h"

using namespace std;

int main(){
	int myrank=4, nb_procs=8, nb_procs_x=2, nx=10, ny=10;
	for(int j=0;j<nb_procs;++j){
		decomposition D(j, nb_procs, nb_procs_x, nx, ny);
		int* index = D.get_index_global_bottom(), length=D.get_myNx();
		cout << "#" << j << ". " << length << ": ";
		for(int i=0;i<length;++i)
			cout << index[i] << " ";
		cout << endl;
	}
	return 0;
}