#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "decomposition.h"

using namespace std;

int main(){
	int myrank=3, nb_procs=4, nb_procs_x=4, nx=17, ny=16;
	decomposition D(myrank, nb_procs, nb_procs_x, nx, ny);
	int* index = D.get_index_x();
	for(int i=0;i<27;++i)
		cout << index[i] << endl;
	return 0;
}