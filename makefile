hi_galaxy: main.o relations.o in_out.o makefile
	g++ -o hi_galaxy -L/usr/local/lib main.o relations.o in_out.o -lgsl -lgslcblas -lm -lCCfits -lcfitsio

main.o: main.cpp relations.cpp in_out.cpp
	g++ -c main.cpp

relations.o: relations.cpp in_out.cpp
	g++ -c relations.cpp

in_out.o: in_out.cpp relations.cpp
	g++ -c in_out.cpp	
