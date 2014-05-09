#include "NS.h"
#include <ZNAC/basic/general.h>

using namespace NS;
using namespace ZNAC;

Diffusion::Diffusion(double Dx, double Dy, unsigned int Nx, unsigned int Ny, unsigned int margin)
	:n(xi, yi, Nx, Ny, margin), Dx(Dx), Dy(Dy)
{
	assert(margin > 0);
}

double &Diffusion::operator()(unsigned int i, unsigned int j)
{
	static double dammy;
	return (dammy = (i == j?-2*(Dx + Dy):ABS((int)(i - j)) == 1?Dx:ABS((int)(i - j)) == n.Nx + 2*n.margin?Dy:0));
}

const double &Diffusion::operator()(unsigned int i, unsigned int j)const
{
	static double dammy;
	return (dammy = (i == j?-2*(Dx + Dy):ABS((int)(i - j)) == 1?Dx:ABS((int)(i - j)) == n.Nx + 2*n.margin?Dy:0));
}

void Diffusion::operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod) const
{
	assert(dom.N() == cod.N() && dom.N() == (n.Nx + 2*n.margin)*(n.Ny + 2*n.margin));

	for(yi = 0; yi < n.Ny; ++yi)
		for(xi = 0; xi < n.Nx; ++xi)
			cod[n] = Dx*(dom[n(1, 0)] + dom[n(-1, 0)]) + Dy*(dom[n(0, 1)] + dom[n(0, -1)]) - 2*(Dx + Dy)*dom[n];

	for(yi = -n.margin; yi < 0; ++yi)
		for(xi = -n.margin; xi < n.Nx + n.margin; ++xi)
			cod[n] = 0;

	for(yi = 0; yi < n.Ny; ++yi)
	{
		for(xi =-n.margin; xi < 0; ++xi)
			cod[n] = 0;
		for(xi = n.Nx; xi < n.Nx + n.margin; ++xi)
			cod[n] = 0;
	}

	for(yi = n.Ny; yi < n.Ny +n.margin; ++yi)
		for(xi = -n.margin; xi < n.Nx + n.margin; ++xi)
			cod[n] = 0;
}

void Advection::operator()(const NS_Grid &u, const NS_Grid &v, NS_Grid &du, NS_Grid &dv)
{
	assert(u.Nx == v.Nx + 1 && u.Ny + 1 == v.Ny);
	assert(u.Nx == du.Nx && u.Ny == du.Ny && v.Nx == dv.Nx && v.Ny == dv.Ny);

	int xi, yi;
	NS_Index un(xi, yi, u.Nx, u.Ny, u.margin), vn(xi, yi, v.Nx, v.Ny, v.margin);

	for(yi = 0; yi < u.Ny; ++yi)
		for(xi = 0; xi < v.Nx; ++xi)
		{
			du[un] = u[un]*dx_inv*(u[un(1, 0)] - u[un(-1, 0)]) + 0.25*((v[vn] + v[vn(-1, 1)]) + (v[vn(0, 1)] + v[vn(-1, 0)]))*dy_inv*(u[un(0, 1)] - u[un(0, -1)]);
			dv[vn] = 0.25*((u[un] + u[un(1, -1)]) + (u[un(1, 0)] + u[un(0, -1)]))*dx_inv*(v[vn(1, 0)] - v[vn(-1, 0)]) + v[vn]*dy_inv*(v[vn(0, 1)] - v[vn(0, -1)]);
		}
}

void Gradient::operator()(const NS_Grid &p, NS_Grid &dpu, NS_Grid &dpv)
{
	assert(p.Nx + 1 == dpu.Nx && p.Ny + 1 && dpv.Ny);
	assert(p.Ny == dpu.Ny && p.Nx == dpv.Nx);

	int xi, yi;
	NS_Index un(xi, yi, dpu.Nx, dpu.Ny, dpu.margin), vn(xi, yi, dpv.Nx, dpv.Ny, dpv.margin), pn(xi, yi, p.Nx, p.Ny, p.margin);
	for(yi = 0; yi < dpv.Ny; ++yi)
		for(xi = 0; xi < dpu.Nx; ++xi)
		{
			dpu[un] = dx_inv*(p[pn] - p[pn(-1, 0)]);
			dpv[vn] = dy_inv*(p[pn] - p[pn(0, -1)]);
		}
}

void Divergence::operator()(const NS_Grid &u, const NS_Grid &v, NS_Grid &p)
{
	assert(p.Nx + 1 == u.Nx && p.Ny == u.Ny);
	assert(p.Nx == v.Nx && p.Ny + 1 == v.Ny);

	int xi, yi;
	NS_Index un(xi, yi, u.Nx, u.Ny, u.margin), vn(xi, yi, v.Nx, v.Ny, v.margin), pn(xi, yi, p.Nx, p.Ny, p.margin);

	for(yi = 0; yi < p.Ny; ++yi)
		for(xi = 0; xi < p.Nx; ++xi)
			p[pn] = dx_inv*(u[un(1, 0)] - u[un]) + dy_inv*(v[vn(0, 1)] - v[vn]);
}

NS_Grid::NS_Grid(unsigned int Nx, unsigned int Ny, unsigned int margin)
	:Nx(Nx), Ny(Ny), margin(margin), buf(new double[(Nx + 2*margin)*(Ny + 2*margin)])
{}

NS_Grid::~NS_Grid(){delete [] buf;}

double &NS_Grid::operator()(int xi, int yi){return buf[(yi + margin)*(Nx + 2*margin) + xi + margin];}
double &NS_Grid::operator[](unsigned int i){assert(i < N()); return buf[i];}
const double &NS_Grid::operator[](unsigned int i) const {assert(i < N()); return buf[i];}
unsigned int NS_Grid::N() const {return (Nx + 2*margin)*(Ny + 2*margin);}

