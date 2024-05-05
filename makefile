#-------------------------------------------------------------------------#
#                         Makefile for pyaneti                            #
#                       Oscar Barragán, Oct. 2016                         #
#-------------------------------------------------------------------------#
#    make       -> compiles the code in its sequential configuration      #
#    make para  -> compiles the code in its parallel configuration        #
#    make clean -> removes the pyaneti.so file                            #
#-------------------------------------------------------------------------#

# let users optionally override f2py, as it could 
# vary, depending on different OSes / environments
ifndef FP
        FP=f2py
endif
#FP=f2py3.8
#FP=f2py2.7
fc=gnu95 
cc=unix 

FLAGS_OMP= -c --quiet --f90flags='-fopenmp' -m
FLAGS = -c --quiet -m
BLIBS = -llapack -lblas
LGOMP = -lgomp

SOURCES=src/constants.f90\
	src/todo.f90\
	src/qpower2.f90\
	src/quad.f90\
        src/ftr.f90\
	src/frv.f90\
	src/bayesian.f90\
	src/matrices.f90\
	src/kernels.f90\
	src/mcmc.f90\
	src/multi-gp-routines.f90\

EXECUTABLE=pyaneti

all: $(EXECUTABLE) 
%.mod: %.90
	$(FC) -c $(SOURCES)
$(EXECUTABLE):$(SOURCES)
	${FP} ${FLAGS} $(EXECUTABLE) $(SOURCES)  --fcompiler=$(fc) ${BLIBS} --compiler=$(cc)


para: $(EXECUTABLE)
	${FP} ${FLAGS_OMP} $(EXECUTABLE) $(SOURCES)  --fcompiler=$(fc) ${BLIBS} ${LGOMP}  --compiler=$(cc)

clean:
	rm *.so src/*.mod src/*.o
