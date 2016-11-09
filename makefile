#-------------------------------------------------------------------------#
#                         Makefile for pyaneti                            #
#                       Oscar Barragán, Oct. 2016                         #
#-------------------------------------------------------------------------#
#    make       -> compiles the code in its sequential configuration      #
#    make para  -> compiles the code in its parallel configuration        #
#    make clean -> removes the pyaneti.so file                            #
#-------------------------------------------------------------------------#

FP=f2py
fc=gnu95 
cc=unix 

FLAGS_OMP= -c -m --f90flags='-fopenmp' 
FLAGS= -c -m 
LGOMP=-L/usr/lib64/ -lgomp
SOURCES= src/ftr.f90\
	src/frv.f90\
	src/quad.f90\
	src/fit_all.f90\
	src/todo.f90

EXECUTABLE=pyaneti

all: $(EXECUTABLE) 
$(EXECUTABLE):$(SOURCES)
	${FP} ${FLAGS} $(EXECUTABLE) $(SOURCES)  --fcompiler=$(fc) --compiler=$(cc)

para: $(EXECUTABLE)
	${FP} ${FLAGS_OMP} $(EXECUTABLE) $(SOURCES)  --fcompiler=$(fc) ${LGOMP}  --compiler=$(cc)

clean:
	rm $(EXECUTABLE).so
