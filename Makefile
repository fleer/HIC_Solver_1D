hydro: HYDRO.cpp
	$(CXX) -std=c++11 -lgsl -lgslcblas -lm -fopenmp HYDRO.cpp
