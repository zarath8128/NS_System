#include "NS.h"

namespace NS
{
	void NS_Gradient(const ZNAC::LA::IVector<double> &p, NS_Vector &cod, double dx)
	{
		const double dx_inv = 2/dx;
		for(unsigned int ny = 1; ny < cod.Ny - 1; ++ny)
			for(unsigned int nx = 1; nx < cod.Nx - 1; ++nx)
			{
				cod.u[ny*cod.Nx + nx] = dx_inv*(p[ny*cod.Nx + nx + 1] - p[ny*cod.Nx + nx - 1]);
				cod.v[ny*cod.Nx + nx] = dx_inv*(p[(ny + 1)*cod.Nx + nx] - p[(ny - 1)*cod.Nx + nx]);
			}
	}

	void NS_Divergence(const NS_Vector &dom, ZNAC::LA::IVector<double> &cod, double dx)
	{
		const double dx_inv = 2/dx;
		for(unsigned int ny = 1; ny < dom.Ny - 1; ++ny)
			for(unsigned int nx = 1; nx < dom.Nx - 1; ++nx)
				cod[ny*dom.Nx + nx] = dx_inv*((dom.u[ny*dom.Nx + nx + 1] - dom.u[ny*dom.Nx - 1]) + (dom.v[(ny + 1)*dom.Nx + nx] - dom.u[(ny - 1)*dom.Nx]));
	}
}
