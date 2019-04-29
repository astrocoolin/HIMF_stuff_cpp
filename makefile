test: NCC.cpp main2.o relations.o makefile NCC.o
	g++ -o test -L/usr/local/lib main2.o -lgsl -lgslcblas -lm -fopenmp

hi_galaxy: main.o relations.o makefile
	g++ -o hi_galaxy -L/usr/local/lib main.o relations.o -lgsl -lgslcblas -lm

main.o: main.cpp relations.cpp 
	g++ -c main.cpp

main2.o: main2.cpp relations.cpp NCC.cpp makefile
	g++ -c main2.cpp -fopenmp

relations.o: relations.cpp NCC.cpp 
	g++ -c relations.cpp

