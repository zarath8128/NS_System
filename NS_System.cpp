#include "NS.h"
#include <iostream>
#include <iomanip>
#include <cstring>

using namespace NS;

NS_System::NS_System(double Re, double a, unsigned int N, double width)
	:Nx(N), Ny(N), Re(Re), a(a), min(-width), max(width), dx((max - min)/(N - 1)), u(Nx, Ny), p(Nx*Ny)
{
	for(unsigned int j = 0; j < Nx; ++j)
		for(unsigned int i = 0; i < Ny; ++i)
			u.u[j*Nx + i] = -(j == Ny - 1?a*((1 - x(i))*(1 + x(i))):0), u.v[j*Nx + i] = 0, p[j*Nx + i] = 0;
}

double NS_System::x(int i){return ((min + i*dx) + (max - (Nx - 1 - i)*dx))*0.5;}
double NS_System::y(int j){return ((min + j*dx) + (max - (Ny - 1 - j)*dx))*0.5;}

NS_Option::NS_Option(int argc, char *argv[])
	:Re(50), a(1), N(64), width(1), arrow_len(2*width/(N*a)),
	Advection(true), Diffusion(true), Force(true), Pressure(true), Auto(true),
	dt((2*width/(N - 1))*(2*width/(N - 1))*0.4), wnd_width(400)
{
	int i = 0;
	if(argc != 1)
		while(++i != argc)
		{
			if(!strcmp(argv[i++], "-Re"))
				if(EOF == sscanf(argv[i], "%d", &Re));
				{
					std::cout << "Invalid input for Re:" << argv[i] << "\n";
				}	

		}
	std::cout << argv[0] << " use deefault parameter" << std::endl;
	print();
}

#define BOOLSTR(b) (b?"true":"false")

constexpr unsigned int walllen = 96;
constexpr unsigned int col1 = 10;
constexpr unsigned int col2 = 10;

template<class T>
void printline(const char *name, const T &val, const char * description)
{
	std::cout << std::setw(col1) << std::setfill(' ') << std::left << name << ":" << std::setw(col2) << std::setfill(' ') << std::right << val << ":\t" << description << "\n";
}

void NS_Option::print()
{
	std::cout << std::setw(walllen) << std::setfill('-') << "\n";
	printline("Re", Re, 				"Reynolds Number.                default = 50");
	printline("a", a, 				"Palaboric Parameter.            default = 1");
	printline("N", N, 				"Number of Dividing space.       default = 64");
	printline("width", width, 			"Length from origin to boundary. default = 1");
	printline("Advection", BOOLSTR(Advection),	"Enable Advection term.          default = true");
	printline("Diffusion", BOOLSTR(Diffusion),	"Enable Diffusion term.          default = true");
	printline("Force", BOOLSTR(Force),		"Enable Force term.              default = true");
	printline("Pressure", BOOLSTR(Pressure),	"Enable Pressure term.           default = true");
	if(!Pressure)
		std::cout << "===============WARNING!! Auto error modifying will NOT operate!!==============" << "\n";
	printline("Auto", BOOLSTR(Auto),		"Enable Auto Progression.        default = true");
	printline("wnd_width", wnd_width,		"Window size(mm).                default = 400");
	std::cout << std::setw(walllen) << std::setfill('-') << "\n";
}
