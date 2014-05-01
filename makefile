ALL=Cavity_glsc test
.POHONY:all, clean

LDLIBS+=-lglsc -lX11
CXXFLAGS=-std=c++0x

all:${ALL}
clean:
	rm -rf ${ALL} *.o

Cavity_glsc:NS.h NS_Vector.o NS_Matrix.o NS_System.o
test:NS.h NS_Vector.o NS_Matrix.o NS_System.o NS_Staggered.o
