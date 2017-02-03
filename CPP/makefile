CC=mpicc

congrad_test.exe: congrad_test.o conjgrad.o comm_ctrl.o decomposition.o operator_matrix.o
	$(CC) -Wall $^ -o $@ -lm

jacobi_test.exe: jacobi_test.o conjgrad.o comm_ctrl.o decomposition.o operator_matrix.o
	$(CC) -Wall $^ -o $@ -lm

	
congrad_test.o: congrad_test.cpp conjgrad.cpp comm_ctrl.cpp decomposition.cpp operator_matrix.cpp congrad_test.h conjgrad.h comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) -Wall -c $@ congrad_test.cpp 
	
jacobi_test.o: jacobi_test.cpp jacobi_test.h jacobi.h comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) -Wall -c $@ jacobi_test.cpp 
	
congrad.o: conjgrad.cpp conjgrad.h comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) -Wall -c $@ conjgrad.cpp 
	
jacobi.o: jacobi.cpp jacobi.h comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) -Wall -c $@ jacobi.cpp 
	
comm_ctrl.o: comm_ctrl.cpp comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) -Wall -c $@ comm_ctrl.cpp
	
decomposition.o: decomposition.cpp decomposition.h 
	$(CC) -Wall -c $@ decomposition.cpp
	
operator_matrix.o: operator_matrix.cpp operator_matrix.h
	$(CC) -Wall -c $@ operator_matrix.cpp

clean: 
	rm -f *.exe *.o
	rm -rf *~