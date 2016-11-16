CC=mpicc
EXEC= main.out

all: $(EXEC) 

main.out: main.o conjgrad.o tools.o cdt_bords.o
	$(CC) -W -Wall $^ -o $@ -lm

%.out: %.o
	$(CC) -W -Wall $< -o $@ -lm

%.o: %.cpp
	$(CC) -W -Wall -c $< -lm

clean: 
	rm -f *.out sol* time* speedup *.o
