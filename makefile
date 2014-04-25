ALL=test
.POHONY:all, clean

LDLIBS+=-lglsc -lX11
CXXFLAGS=-std=c++0x

all:${ALL}
clean:
	rm -rf ${ALL}

test:NS.h NS_Vector.cpp NS_Matrix.cpp NS_Function.cpp NS_System.cpp
