FC = gfortran
LINK = $(FC)
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

all: mat_mult yeval demo2 fine_grain coarse_grain

mat_mult: mod_rand.o matrix_multiply.o 
	$(LINK) $(LFLAGS) $^ -o mat_mult

yeval: yeval.o
	$(LINK) $(LFLAGS) $^ -o yeval

demo2: demo2.o
	$(LINK) $(LFLAGS) $^ -o demo2

fine_grain: fine_grain.o
	$(LINK) $(LFLAGS) $^ -o fine_grain

coarse_grain: coarse_grain.o
	$(LINK) $(LFLAGS) $^ -o coarse_grain

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

clean:
	-rm -f mod_rand.o mod_rand.mod matrix_multiply.o mat_mult
	-rm -f yeval.o yeval
	-rm -f demo2.o demo2
	-rm -f fine_grain.o fine_grain
	-rm -f coarse_grain.o coarse_grain

new:
	$(MAKE) clean
	$(MAKE) all


### DO NOT remove this line - make depends on it ###