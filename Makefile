hydro: HYDRO.cpp
	$(CXX) -o hydro -std=c++11 -lgsl -lgslcblas -lm -fopenmp HYDRO.cpp
