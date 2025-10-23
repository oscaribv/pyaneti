#-------------------------------------------------------------------------#
#                         Makefile for pyaneti                            #
#                       Oscar Barragán, Oct. 2016                         #
#-------------------------------------------------------------------------#
#    make                   -> compiles the code with GNU compilers       #
#    make COMPILER=intel    -> compiles the code with Intel compilers     #
#    make para              -> parallel build with GNU compilers          #
#    make para COMPILER=intel -> parallel build with Intel compilers  #
#    make clean             -> removes the pyaneti.so file                #
#-------------------------------------------------------------------------#

# Default compiler is GNU
COMPILER ?= gnu

# Compiler settings
ifeq ($(COMPILER),gnu)
    FP = python -m numpy.f2py
    FC = gfortran
    CC = unix
    FLAGS_OMP = -c --quiet --f90flags='-fopenmp' -m
    FLAGS = -c --quiet -m
    BLIBS = -llapack -lblas
    LGOMP = -lgomp
else ifeq ($(COMPILER),intel)
    FP = f2py3.11
    FC = ifx
    CC = icx-cc
    FLAGS_OMP = -c -m --quiet --f90flags='-qopenmp'
    FLAGS = -c -m
    BLIBS = -mkl
    LGOMP = -qopenmp
else
    $(error Unknown compiler specified. Use COMPILER=gnu or COMPILER=intel)
endif

# Source files
SOURCES = src/constants.f90 \
          src/todo.f90 \
          src/qpower2.f90 \
          src/quad.f90 \
          src/ftr.f90 \
          src/frv.f90 \
          src/bayesian.f90 \
          src/matrices.f90 \
          src/kernels.f90 \
          src/mcmc.f90 \
          src/multi-gp-routines.f90

# Output binary
EXECUTABLE = pyaneti

# Default rule (sequential build)
all: $(EXECUTABLE)

# Object files
OBJS = $(SOURCES:.f90=.o)

# Dependency files
DEPS = $(SOURCES:.f90=.d)

# Compilation rule
$(EXECUTABLE): $(OBJS)
	$(FP) $(FLAGS) $(EXECUTABLE) $(SOURCES) --fcompiler=$(FC) $(BLIBS) --compiler=$(CC)

# Parallel compilation rule
para: $(OBJS)
	$(FP) $(FLAGS_OMP) $(EXECUTABLE) $(SOURCES) --fcompiler=$(FC) $(BLIBS) $(LGOMP) --compiler=$(CC)

# Compile source files to object files
%.o: %.f90
	$(FC) -c $< -o $@

# Include dependencies
-include $(DEPS)

# Clean rule
clean:
	rm -f $(EXECUTABLE)*so src/*.mod src/*.o

# Phony targets
.PHONY: all parallel clean

