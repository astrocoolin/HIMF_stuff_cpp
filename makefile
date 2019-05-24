observe: NCC.h main.o  makefile 
	g++ -o observe -L/usr/local/lib main.o -lgsl -lgslcblas -lm -fopenmp

main.o: main.cpp relations.h NCC.h makefile
	g++ -Wall -O3 -c main.cpp -fopenmp
