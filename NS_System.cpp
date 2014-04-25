#include "NS.h"

using namespace NS;

NS_System::NS_System(double Re, double a, unsigned int N, double width)
	:Nx(N), Ny(N), Re(Re), a(a), min(-width), max(width), dx((max - min)/(N - 1)), u(Nx, Ny), p(Nx*Ny)
{
	for(unsigned int j = 0; j < N; ++j)
		for(unsigned int i = 0; i < N; ++i)
			u.u[j*Nx + i] = (j == 0?a*(1 - x(i))*(1 + x(i)):0), u.v[j*Nx + i] = 0, p[j*Nx + i] = 0;
}

double NS_System::x(int i){return min + i*dx;}
double NS_System::y(int j){return max - j*dx;}
