hi_galaxy: main.o relations.o makefile
	g++ -o hi_galaxy -L/usr/local/lib main.o -lgsl -lgslcblas -lm

main.o: main.cpp relations.cpp
	g++ -c main.cpp

relations.o: relations.cpp
	g++ -c relations.cpp
