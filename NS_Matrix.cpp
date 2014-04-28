#include "NS.h"
#include <ZNAC/basic/general.h>

using namespace NS;
using namespace ZNAC;

NS_Matrix::NS_Matrix(unsigned int Nx, unsigned int Ny):Nx(Nx), Ny(Ny){}

NS_Diffusion::NS_Diffusion(unsigned int Nx, unsigned int Ny, double Re, double dx):NS_Matrix(Nx, Ny), Redx2_inv(1/(Re*dx*dx)){}
double &NS_Diffusion::operator()(unsigned int r, unsigned int c)
{
	static double dammy;
	return r == c?dammy = -4:ABS(r - c) == 1 || ABS(r - c) == Nx?dammy = 1:dammy = 0;
}
const double &NS_Diffusion::operator()(unsigned int r, unsigned int c)const
{
	static double dammy;
	return r == c?dammy = -4:ABS(r - c) == 1 || ABS(r - c) == Nx?dammy = 1:dammy = 0;
}

void NS_Diffusion::operator()(const LA::IVector<double> &dom, LA::IVector<double> &cod) const
{
	for(unsigned int ny = 1; ny < Ny - 1; ++ny)
		for(unsigned int nx = 1; nx < Nx - 1; ++nx)
		{
			cod[ny*Nx +nx] = Redx2_inv*(-4*dom[ny*Nx + nx] + ((dom[(ny + 1)*Nx + nx] + dom[(ny - 1)*Nx + nx]) + (dom[ny*Nx + nx + 1] + dom[ny*Nx + nx - 1])));
			cod[ny*Nx +nx + Nx*Ny] = Redx2_inv*(-4*dom[ny*Nx + nx + Nx*Ny] + ((dom[(ny + 1)*Nx + nx + Nx*Ny] + dom[(ny - 1)*Nx + nx + Nx*Ny]) + (dom[ny*Nx + nx + 1 + Nx*Ny] + dom[ny*Nx + nx - 1 + Nx*Ny])));
		}
}

NS_Diffusion_S::NS_Diffusion_S(unsigned int Nx, unsigned int Ny, double dx):NS_Matrix(Nx, Ny), dx2_inv(1/(dx*dx)){}
double &NS_Diffusion_S::operator()(unsigned int r, unsigned int c)
{
	static double dammy;
	return r == c?dammy = 4:ABS(r - c) == -1 || ABS(r - c) == Nx?dammy = -1:dammy = 0;
}
const double &NS_Diffusion_S::operator()(unsigned int r, unsigned int c)const
{
	static double dammy;
	return r == c?dammy = 4:ABS(r - c) == -1 || ABS(r - c) == Nx?dammy = -1:dammy = 0;
}

void NS_Diffusion_S::operator()(const LA::IVector<double> &dom, LA::IVector<double> &cod) const
{
//	for(unsigned int n = 0; n < Nx; ++n)
//	{
//		cod[n] = cod[(Ny - 1)*Nx + n] = cod[n*Nx] = cod[n*Nx + Nx - 1] = 0;
//	}
	for(unsigned int ny = 1; ny < Ny - 1; ++ny)
		for(unsigned int nx = 1; nx < Nx - 1; ++nx)
			cod[ny*Nx +nx] = dx2_inv*(4*dom[ny*Nx + nx] - ((dom[(ny + 1)*Nx + nx] + dom[(ny - 1)*Nx + nx]) + (dom[ny*Nx + nx + 1] + dom[ny*Nx + nx - 1])));
}

NS_Advection::NS_Advection(unsigned int Nx, unsigned int Ny, double dx):NS_Matrix(Nx, Ny), dx_inv(1/(2*dx)){}
double &NS_Advection::operator()(unsigned int r, unsigned int c)
{
	static double dammy;
	return (r - c) == 1 || (r - c) == Nx?dammy = 1:r - c == -1 || r - c == -Nx?dammy = -1:dammy = 0;
}
const double &NS_Advection::operator()(unsigned int r, unsigned int c)const
{
	static double dammy;
	return (r - c) == 1 || (r - c) == Nx?dammy = 1:r - c == -1 || r - c == -Nx?dammy = -1:dammy = 0;
}

void NS_Advection::operator()(const LA::IVector<double> &dom, LA::IVector<double> &cod) const
{
	for(unsigned int ny = 1; ny < Ny - 1; ++ny)
		for(unsigned int nx = 1; nx < Nx - 1; ++nx)
		{
			cod[ny*Nx +nx] = dx_inv*(dom[ny*Nx + nx + Nx*Ny]*(dom[(ny + 1)*Nx + nx] - dom[(ny - 1)*Nx + nx]) + dom[ny*Nx + nx]*(dom[ny*Nx + nx + 1] - dom[ny*Nx + nx - 1]));
			cod[ny*Nx +nx + Nx*Ny] = dx_inv*(dom[ny*Nx + nx + Nx*Ny]*(dom[(ny + 1)*Nx + nx + Nx*Ny] - dom[(ny - 1)*Nx + nx + Nx*Ny]) + dom[ny*Nx + nx]*(dom[ny*Nx + nx + 1 + Nx*Ny] - dom[ny*Nx + nx - 1 + Nx*Ny]));
		}
}
