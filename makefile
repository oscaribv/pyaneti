FP=f2py
fc=gnu95 
cc=unix 

FLAGS= -c -m --f90flags='-fopenmp' 
#FLAGS= -c -m 
LIB=-L/usr/lib64/ -lgomp
SOURCES= src/ftr.f90\
	src/frv.f90\
	src/quad.f90\
	src/todo.f90

EXECUTABLE=pyaneti

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	${FP} ${FLAGS} $(EXECUTABLE) $(SOURCES)  --fcompiler=$(fc) ${LIB}  --compiler=$(cc)

clean:
	rm $(EXECUTABLE).so
