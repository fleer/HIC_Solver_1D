hydro: HYDRO.cpp
	$(CXX) -std=c11 -lgsl -lgslcblas -lm -fopenmp HYDRO.cpp
