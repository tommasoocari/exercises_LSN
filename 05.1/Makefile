CC = g++
CFLAGS = -Wall -O3 --std=c++11
PRNG = ../Parallel_Random_Number_Generator/
LIB = ../lib/

main.exe : main.o $(PRNG)random.o $(LIB)error.o $(LIB)lib.o distribution.o
	$(CC) $(PRNG)random.o main.o $(LIB)error.o $(LIB)lib.o distribution.o -o main.exe
main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : $(PRNG)random.cpp $(PRNG)random.h
	$(CC) -c $(PRNG)random.cpp -o $(PRNG)random.o $(CFLAGS)
error.o : $(LIB)error.cpp $(LIB)error.h
	$(CC) -c $(LIB)error.cpp -o $(LIB)error.o $(CFLAGS)
lib.o : $(LIB)lib.cpp $(LIB)lib.h
	$(CC) -c $(LIB)lib.cpp -o $(LIB)lib.o $(CFLAGS)
distribution.o : distribution.cpp distribution.h
	$(CC) -c distribution.cpp -o distribution.o $(CFLAGS)
clean :
	rm *.o main.exe
