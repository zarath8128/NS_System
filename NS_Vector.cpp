#include "NS.h"
#include <iostream>

using namespace NS;

NS_Vector::NS_Vector(unsigned int Nx, unsigned int Ny):buf(new double[2*Nx*Ny]), Nx(Nx), Ny(Ny), u(buf), v(buf + Nx*Ny){}
NS_Vector::~NS_Vector(){delete [] buf;}

double &NS_Vector::operator[](unsigned int i){return buf[i];}
const double &NS_Vector::operator[](unsigned int i) const {return buf[i];}

unsigned int NS_Vector::N() const {return Nx*Ny*2;}
