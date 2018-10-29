hi_galaxy: main.o relations.o
	g++ -o hi_galaxy main.o

main.o: main.cpp relations.cpp
	g++ -c main.cpp

relations.o: relations.cpp
	g++ -c relations.cpp
