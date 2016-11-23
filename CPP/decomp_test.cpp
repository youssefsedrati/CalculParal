#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "decomposition.h"

using namespace std;

int main(){
	int myrank=1, nb_procs=3, nb_procs_x=3, nx=10, ny=10;
	decomposition D(myrank, nb_procs, nb_procs_x, nx, ny);
	int* index = D.get_index_global();
	for(int i=0;i<10;++i)
		cout << index[i] << endl;
	return 0;
}