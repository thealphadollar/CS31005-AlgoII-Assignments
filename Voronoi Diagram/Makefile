a.out: main.o VoronoiDiagram.o
	g++ main.o VoronoiDiagram.o
main.o: main.cpp VoronoiDiagram.o
	g++ -c main.cpp
VoronoiDiagram.o: VoronoiDiagram.cpp VoronoiDiagram.hpp
	g++ -c 	VoronoiDiagram.cpp
clean:
	rm a.out main.o VoronoiDiagram.o t4.svg
