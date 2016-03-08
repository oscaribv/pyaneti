FP=f2py
fc=gfortran
cc=unix 

FLAGS= -c -m
SOURCES= 	ftr.f90\
					frv.f90\
					quad.f90\
					todo.f90

EXECUTABLE=pyaneti

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCES)
	${FP} ${FLAGS} $(EXECUTABLE) $(SOURCES) --fcompiler=$(fc) --compiler=$(cc)

clean:
	rm $(EXECUTABLE).so
