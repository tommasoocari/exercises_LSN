CC = g++
CFLAGS = -Wall -O3 --std=c++11
PRNG = ../Parallel_Random_Number_Generator/

TSP.exe : TSP.o $(PRNG)random.o
	$(CC) $(PRNG)random.o TSP.o -o TSP.exe
TSP.o : TSP.cpp TSP.h
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)
random.o : $(PRNG)random.cpp $(PRNG)random.h
	$(CC) -c $(PRNG)random.cpp -o $(PRNG)random.o $(CFLAGS)
clean :
	rm *.o TSP.exe
