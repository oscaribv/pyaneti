FP=f2py
fc=gnu95 
cc=unix 

FLAGS_OMP= -c -m --f90flags='-fopenmp -lgomp' 
FLAGS= -c -m 
SOURCES= src/ftr.f90\
	src/frv.f90\
	src/quad.f90\
	src/todo.f90

EXECUTABLE=pyaneti

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	${FP} ${FLAGS} $(EXECUTABLE) $(SOURCES) --fcompiler=$(fc) --compiler=$(cc)

para: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	${FP} ${FLAGS_OMP} $(EXECUTABLE) $(SOURCES) --fcompiler=$(fc) --compiler=$(cc)

clean:
	rm $(EXECUTABLE).so
