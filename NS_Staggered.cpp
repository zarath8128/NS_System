#include "NS.h"

using namespace NS;

NS_Staggered::Coor::Coor(double dx, double offset, int &i):dx(dx), offset(offset), i(i){}
NS_Staggered::Coor::operator double()const{return offset + dx*i;}

NS_Staggered::Index::Index(int &i, int &j, unsigned int Nx, unsigned int Ny, unsigned int margin):i(i), j(j), Nx(Nx), Ny(Ny), margin(margin){}
unsigned int NS_Staggered::Index::operator()(int i, int j){return (j + margin)*(Nx + margin*2) + i + margin;}
NS_Staggered::Index::operator unsigned int(){return (Nx + margin*2)*(j + margin) + i + margin;}

NS_Staggered::NS_Staggered(double Re, unsigned int Nx, unsigned int Ny,  
		double (*u)(double x, double y, unsigned int i, unsigned int j), 
		double (*v)(double x, double y, unsigned int i, unsigned int j), 
		double (*p)(double x, double y, unsigned int i, unsigned int j), 
	       	double left, double right, double bottom, double top, unsigned int margin)
	:Nx(Nx), Ny(Ny), margin(margin), Re(Re), left(left), right(right), bottom(bottom), top(top), dx((right - left)/Nx), dy((top - bottom)/Ny), 
	redx2_inv(1./(Re*dx*dx)), redy2_inv(1./(Re*dy*dy)), dx_inv(1/dx), dy_inv(1/dy),
	u(Nx + 1 + 2*margin, Ny + 2*margin), p((Nx + 2*margin)*(Ny + 2*margin)), xoffset(dx*.5), yoffset(dy*.5), 
	i(0), j(0), un(i, j, Nx + 1, Ny, margin), vn(i, j, Nx, Ny + 1, margin), pn(i, j, Nx, Ny, margin),
	ux(dx, left, i), uy(dy, yoffset + bottom, j), vx(dx, xoffset + left, i), vy(dy, bottom, j), px(dx, xoffset + left, i), py(dy, yoffset + bottom, j),
	diffusion_p(Nx, Ny, margin, dx, dy)
{
	for(j = 0; j < Ny + 1; ++j)
		for(i = 0; i < Nx + 1; ++i)
			this->u.u[un] = u(ux, uy, i, j), this->u.v[un] = v(vx, vy, i, j);
	for(j = 0; j < Ny; ++j)
		for(i = 0; i < Nx; ++i)
			this->p[pn] = p(px, py, i, j);
}

void NS_Staggered::Diffusion_u(NS_Vector &du)
{
	for(j = 0; j < Ny; ++j)
		for(i = 0; i < Nx + 1; ++i)
			du.u[un] = redx2_inv*(u.u[un(i + 1, j)] + u.u[un(i - 1, j)] - 2*u.u[un]) + redy2_inv*(u.u[un(i, j - 1)] + u.u[un(i, j + 1)] - 2*u.u[un]);
	for(j = 0; j < Ny + 1; ++j)
		for(i = 0; i < Nx; ++i)
			du.v[vn] = redx2_inv*(u.v[vn(i + 1, j)] + u.v[vn(i - 1, j)] - 2*u.v[vn]) + redy2_inv*(u.v[vn(i, j - 1)] + u.v[vn(i, j + 1)] - 2*u.v[vn]);
}


void NS_Staggered::Gradient(const ZNAC::LA::IVector<double> &s, NS_Vector &ds)
{
	for(j = 0; j < Ny; ++j)
		for(i = 0; i < Nx + 1; ++i)
			ds.u[un] = dx_inv*(s[pn] - s[pn(i - 1, j)]);
	for(j = 0; j < Ny + 1; ++j)
		for(i = 0; i < Nx; ++i)
			ds.v[vn] = dy_inv*(s[pn] - s[pn(i, j - 1)]);
}

void NS_Staggered::Advection(NS_Vector &du)
{
	for(j = 0; j < Ny; ++j)
		for(i = 0; i < Nx + 1; ++i)
			du.u[un] = u.u[un]*dx_inv*(u.u[un(i + 1, j)] - u.u[un(i - 1, j)]) + u.v[vn]*dy_inv*(u.u[un(i, j - 1)] - u.u[un(i, j + 1)]);
	for(j = 0; j < Ny + 1; ++j)
		for(i = 0; i < Nx; ++i)
			du.v[vn] = u.u[un]*dx_inv*(u.v[vn(i + 1, j)] - u.v[vn(i - 1, j)]) + u.v[vn]*dy_inv*(u.v[vn(i, j - 1)] - u.v[vn(i, j + 1)]);
}

#include <iostream>

void NS_Staggered::Divergence(ZNAC::LA::IVector<double> &du)
{
	for(j = 0; j < Ny; ++j)
		for(i = 0; i < Nx; ++i)
			du[pn] = dx_inv*(u.u[un(i + 1, j)] - u.u[un]) + dy_inv*(u.v[vn(i, j + 1)] - u.v[vn]);
		
}

NS_Staggered::Diffusion::Diffusion(unsigned int Nx, unsigned int Ny, unsigned int margin, double dx, double dy)
	:Nx(Nx + 2*margin), Ny(Ny + 2*margin), margin(margin), dx2_inv(1./(dx*dx)), dy2_inv(1./(dy*dy)){}
double &NS_Staggered::Diffusion::operator()(unsigned int i, unsigned int j){static double dammy = 0; return dammy = 0;}
const double &NS_Staggered::Diffusion::operator()(unsigned int i, unsigned int j)const {static double dammy = 0; return dammy = 0;}
void NS_Staggered::Diffusion::operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const
{
	for(unsigned int j = 0; j < Ny; ++j)
		for(unsigned int i = 0; i < Nx; ++i)
			cod[j*Nx + i] = (j>=margin && j < Ny - margin && i >= margin && i < Nx - margin)? -dx2_inv*(dom[j*Nx + i - 1] + dom[j*Nx + i + 1] - 2*dom[j*Nx + i]) - dy2_inv*(dom[(j - 1)*Nx + i] + dom[(j + 1)*Nx + i] - 2*dom[j*Nx + i]):0;
}
