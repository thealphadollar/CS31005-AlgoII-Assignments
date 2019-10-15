a.out: main.o Point2D.o VoronoiDiagram.o Beachline.o
	g++ main.o VoronoiDiagram.o
main.o: main.cpp Point2D.o VoronoiDiagram.o Beachline.o
	g++ -c main.cpp
Point2D.o: Point2D.cpp Point2D.h
	g++ -c Point2D.cpp
VoronoiDiagram.o: VoronoiDiagram.cpp VoronoiDiagram.hpp Beachline.o Parabola.o Circle.o DCEL.o Point2D.o
	g++ -c 	VoronoiDiagram.cpp
Beachline.o: Beachline.cpp Beachline.hpp Parabola.o DCEL.o
	g++ -c Beachline.cpp
Parabola.o: Parabola.cpp Parabola.hpp Point2D.o
	g++ -c Parabola.cpp 
Circle.o: Circle.cpp Circle.hpp Point2D.o
	g++ -c Circle.cpp
DCEL.o: DCEL.cpp DCEL.hpp Point2D.o
	g++ -c DCEL.cpp
clean:
	rm a.out main.o Point2D.o VoronoiDiagram.o Beachline.o Parabola.o Circle.o DCEL.o t4.svg
