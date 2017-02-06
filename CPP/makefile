CC=g++
FLAGS = -Wall -c

all: jacobi_test.exe conjgrad_test.exe

jacobi_test.exe: jacobi_test.o jacobi.o conjgrad.o comm_ctrl.o decomposition.o operator_matrix.o
	$(CC) $^ -o $@ -lmsmpi
	
conjgrad_test.exe: conjgrad_test.o conjgrad.o comm_ctrl.o decomposition.o operator_matrix.o
	$(CC) $^ -o $@ -lmsmpi

conjgrad_test.o: conjgrad_test.cpp conjgrad.h comm_ctrl.h decomposition.h operator_matrix.h RHS.h
	$(CC) $(FLAGS) conjgrad_test.cpp 
	
jacobi_test.o: jacobi_test.cpp jacobi.h comm_ctrl.h decomposition.h operator_matrix.h RHS.h
	$(CC) $(FLAGS) jacobi_test.cpp 
	
congrad.o: conjgrad.cpp conjgrad.h comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) $(FLAGS) conjgrad.cpp 
	
jacobi.o: jacobi.cpp jacobi.h comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) $(FLAGS) jacobi.cpp 
	
comm_ctrl.o: comm_ctrl.cpp comm_ctrl.h decomposition.h operator_matrix.h
	$(CC) $(FLAGS) comm_ctrl.cpp
	
decomposition.o: decomposition.cpp decomposition.h 
	$(CC) $(FLAGS) decomposition.cpp
	
operator_matrix.o: operator_matrix.cpp operator_matrix.h
	$(CC) $(FLAGS) operator_matrix.cpp

clean: 
	rm -f *.exe *.o
	rm -rf *~