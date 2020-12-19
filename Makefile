all: omp

omp : coolturkey.cpp
	@mpicxx -o ct.out coolturkey.cpp -fopenmp
	@sbatch sprob1_openmp.sh

clear: 
	rm ct.out
