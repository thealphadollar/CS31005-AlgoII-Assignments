# Assignment 4: Voronoi Diagrams
# Shivam Kumar Jha
# 17CS30033
# Makefile for easy compilation

a.out: main.o Voronoi_Structure.o
	g++ main.o Voronoi_Structure.o
main.o: main.cpp Voronoi_Structure.o
	g++ -c main.cpp
Voronoi_Structure.o: Voronoi_Structure.cpp Voronoi_Structure.hpp
	g++ -c 	Voronoi_Structure.cpp
clean:
	rm a.out main.o Voronoi_Structure.o t4.svg
