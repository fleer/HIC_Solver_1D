all:
	g++ hydro.cpp -o hydro -lgsl -lgslcblas -lm -fopenmp 
	#g++-7 hydro.cpp -o hydro -lgsl -lgslcblas -lm -fopenmp 
	#clang++ hydro.cpp -o hydro -lgsl -lgslcblas -lm -fopenmp 
