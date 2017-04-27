FC = gfortran
MPI_FC = mpif90
LINK = $(FC)
MPI_LINK = $(MPI_FC)
# GFortran flags
FFLAGS = -O3 -funroll-loops -finline-functions -fdefault-real-8 -fopenmp
# Intel Flags
# FFLAGS = -O3 -fast -ipo -no-prec-div -openmp
NAME = mat_mult
TEST_ITER = 1000

# Uncomment the appropriate flags for your platform
# Linux
LFLAGS = $(FFLAGS) -lblas
# Mac OS X framework specific, flag above also works for Mac OS X
# LFLAGS = $(FFLAGS) -framework veclib

# Default fortran rules
%.o : %.f90 ; $(FC) -c $< -o $@ $(FFLAGS)

all: mat_mult yeval demo2 fine_grain coarse_grain hello_world_omp \
	hello_world_mpi compute_pi note_passing \
	jacobi_omp1 jacobi_omp2 jacobi_mpi

# Beginning Demo
mat_mult: mod_rand.o matrix_multiply.o 
	$(LINK) $(LFLAGS) $^ -o $@

mat_mult_test: mat_mult
	-echo "Non-threaded Tests:"
	./mat_mult $(TEST_ITER) 1 1
	./mat_mult $(TEST_ITER) 2 1
	./mat_mult $(TEST_ITER) 3 1
	./mat_mult $(TEST_ITER) 4 1
	-echo "Threaded Tests:"
	./mat_mult $(TEST_ITER) 1 4
	./mat_mult $(TEST_ITER) 2 4
	./mat_mult $(TEST_ITER) 3 4
	./mat_mult $(TEST_ITER) 4 4

# OpenMP Demos
hello_world_omp: hello_world_omp.o
	$(LINK) $(LFLAGS) $^ -o $@

yeval: yeval.o
	$(LINK) $(LFLAGS) $^ -o $@

demo2: demo2.o
	$(LINK) $(LFLAGS) $^ -o $@

fine_grain: fine_grain.o
	$(LINK) $(LFLAGS) $^ -o $@

coarse_grain: coarse_grain.o
	$(LINK) $(LFLAGS) $^ -o $@

jacobi_omp1: jacobi_omp1.o
	$(LINK) $(LFLAGS) $^ -o $@

jacobi_omp2: jacobi_omp2.o
	$(LINK) $(LFLAGS) $^ -o $@

jacobi_omp3: jacobi_omp3.o
	$(LINK) $(LFLAGS) $^ -o $@

# MPI demos
hello_world_mpi.o: hello_world_mpi.f90
	$(MPI_FC) -c $< -o $@ $(FFLAGS)

hello_world_mpi: hello_world_mpi.o
	$(MPI_LINK) $^ -o $@

compute_pi.o: compute_pi.f90
	$(MPI_FC) -c $< -o $@ $(FFLAGS)

compute_pi: compute_pi.o
	$(MPI_LINK) $^ -o $@

note_passing.o: note_passing.f90
	$(MPI_FC) -c $< -o $@ $(FFLAGS)

note_passing: note_passing.o
	$(MPI_LINK) $^ -o $@

jacobi_mpi.o: jacobi_mpi.f90
	$(MPI_FC) -c $< -o $@ $(FFLAGS)

jacobi_mpi: jacobi_mpi.o
	$(MPI_LINK) $^ -o $@

clean:
	-rm -f mod_rand.o mod_rand.mod matrix_multiply.o mat_mult
	-rm -f yeval.o yeval
	-rm -f demo2.o demo2
	-rm -f fine_grain.o fine_grain
	-rm -f coarse_grain.o coarse_grain
	-rm -f hello_world_omp.o hello_world_omp
	-rm -f hello_world_mpi.o hello_world_mpi
	-rm -f compute_pi.o compute_pi
	-rm -f note_passing.o note_passing
	-rm -f jacobi_omp1.o jacobi_omp1 jacobi_omp1.txt
	-rm -f jacobi_omp2.o jacobi_omp2 jacobi_omp2.txt
	-rm -f jacobi_mpi.o jacobi_mpi jacobi_mpi.txt

new:
	$(MAKE) clean
	$(MAKE) all


### DO NOT remove this line - make depends on it ###