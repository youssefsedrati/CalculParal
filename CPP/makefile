CC=mpicc

congrad_test.exe: congrad_test.o conjgrad.o comm_ctrl.o decomposition.o operator_matrix.o
	$(CC) -W -Wall $^ -o $@ -lm

jacobi_test.exe: jacobi_test.o conjgrad.o comm_ctrl.o decomposition.o operator_matrix.o
	$(CC) -W -Wall $^ -o $@ -lm
	
%.o: %.cpp
	$(CC) -W -Wall -c $< -lm

clean: 
	rm -f *.exe *.o
	rm -rf *~