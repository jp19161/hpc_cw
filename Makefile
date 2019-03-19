default: all

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

Model.o: Model.cpp Model.h
	mpicxx -std=c++11 -Wall -o Model.o -c Model.cpp

Burgers.o: Burgers.cpp Burgers.h Model.h
	mpicxx -std=c++11 -Wall -o Burgers.o -c Burgers.cpp

compile: main.o  Model.o Burgers.o
	mpicxx -o HPC_Prog main.o  Model.o Burgers.o -O3 -ffast-math -funroll-loops  -march=native -ftree-vectorize

.PHONY: clean 
	target 

diff: compile
	mpiexec -np 1 HPC_Prog 0 0 0 1 10.0 21 21 4000 1.0 1 1   

advx: compile
	mpiexec -np 1 HPC_Prog 1 0 0 0 10.0 21 21 4000 1.0 1 1  

advy: compile
	mpiexec -np 1 HPC_Prog 0 1 0 0 10.0 21 21 4000 1.0 1 1  

burg: compile
	mpiexec -np 1 HPC_Prog 1.0 0.5 1.0 0.02 10.0 2001 2001 4000 1.0 1 1

diffp: compile
	mpiexec -np 16 HPC_Prog 0 0 0 1 10.0 201 201 1600 1.0 4 4   

advxp: compile
	mpiexec -np 16 HPC_Prog 1 0 0 0 10.0 21 21 4000 1.0 4 4 

advyp: compile
	mpiexec -np 16 HPC_Prog 0 1 0 0 10.0 21 21 4000 1.0 4 4  

burgp: compile
	mpiexec -np 16 HPC_Prog 1.0 0.5 1.0 0.02 10.0 2001 2001 4000 1.0 4 4        

clean:
	-rm -f *.o HPC_Prog   # Clean up (and ignore any errors)

all: diff advx advy burg clean

# ax ay b c L Nx Ny Nt T px py