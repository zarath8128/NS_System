ALL=Cavity_glsc
.POHONY:all, clean

LDLIBS+=-lglsc -lX11
CXXFLAGS=-std=c++0x

all:${ALL}
clean:
	rm -rf ${ALL}

Cavity_glsc:NS.h NS_Vector.cpp NS_Matrix.cpp NS_System.cpp
