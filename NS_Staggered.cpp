#include "NS.h"

using namespace NS;

NS_Staggered::Coor::Coor(double dx, double offset, unsigned int &i):dx(dx), offset(offset), i(i){}
NS_Staggered::Coor::operator double()const{return offset + dx*i;}

NS_Staggered::Index::Index(unsigned int &i, unsigned int &j, unsigned int Nx, unsigned int Ny):i(i), j(j), Nx(Nx), Ny(Ny){}
NS_Staggered::Index::operator unsigned int(){return Nx*j + i;}

NS_Staggered::NS_Staggered(double Re, unsigned int Nx, unsigned int Ny,  double (*u)(double x, double y), double (*v)(double x, double y), double(*p)(double x, double y), double rx, double ry)
	:Nx(Nx), Ny(Ny), Re(Re), xmax(rx), xmin(-rx), ymax(ry), ymin(-ry), dx((xmax - xmin)/Nx), dy((ymax - ymin)/Ny), u(Nx + 1, Ny + 1), p(Nx*Ny), xoffset(dx*.5), yoffset(dy*.5), 
	i(0), j(0), un(i, j, Nx + 1, Ny + 1), pn(i, j, Nx, Ny),
	ux(dx, xmin, i), uy(dy, yoffset + ymin, j), vx(dx, xoffset + xmin, i), vy(dy, ymin, j), px(dx, xoffset + xmin, i), py(dy, yoffset + ymin, j)
{
	for(j = 0; j < Ny + 1; ++j)
		for(i = 0; i < Nx + 1; ++i)
			this->u.u[un] = u(ux, uy), this->u.v[un] = v(vx, vy);
	for(j = 0; j < Ny; ++j)
		for(i = 0; i < Nx; ++i)
			this->p[pn] = p(px, py);
}
