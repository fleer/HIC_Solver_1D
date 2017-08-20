all:
	g++ hydro.cpp -std=c++11 -Wall -o hydro -lgsl -lgslcblas -lm -fopenmp 
