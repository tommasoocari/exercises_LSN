CC = g++
CFLAGS = -Wall -O3 --std=c++11
PRNG = ../Parallel_Random_Number_Generator/

Variational_MC.exe : Variational_MC.o $(PRNG)random.o Wavefunction.o
	$(CC) $(PRNG)random.o Variational_MC.o Wavefunction.o -o Variational_MC.exe
Variational_MC.o : Variational_MC.cpp
	$(CC) -c Variational_MC.cpp -o Variational_MC.o $(CFLAGS)
random.o : $(PRNG)random.cpp $(PRNG)random.h
	$(CC) -c $(PRNG)random.cpp -o $(PRNG)random.o $(CFLAGS)
Wavefunction.o : Wavefunction.cpp Wavefunction.h
	$(CC) -c Wavefunction.cpp -o Wavefunction.o $(CFLAGS)
clean :
	rm *.o Variational_MC.exe
