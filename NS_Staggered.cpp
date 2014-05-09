#include "NS.h"
#include <iostream>

using namespace NS;

NS_Staggered::NS_Staggered(double Re, unsigned int Nx, unsigned int Ny, 
		double (*U)(double x, double y, int i, int j), double (*V)(double x, double y, int i, int j), double (*P)(double x, double y, int i, int j), 
		double left, double right, double bottom, double top, unsigned int margin)
	:Re(Re), dx((right - left)/Nx), dy((top - bottom)/Ny), 
	u(Nx + 1, Ny, margin), v(Nx, Ny + 1, margin), p(Nx, Ny, margin),
	ux(dx, left), uy(dy, bottom + dy*0.5), vx(dx, left + dx*0.5), vy(dy, bottom), px(dx, left + 0.5*dx), py(dy, bottom + 0.5*dy)
{
	for(uy.i = -(int)margin; uy.i < (int)(u.Ny + margin); ++uy.i)
		for(ux.i = -(int)margin; ux.i < (int)(u.Nx + margin); ++ux.i)
			u(ux.i, uy.i) = U(ux, uy, ux.i, uy.i);

	for(vy.i = -(int)margin; vy.i < (int)(v.Ny + margin); ++vy.i)
		for(vx.i = -(int)margin; vx.i < (int)(v.Nx + margin); ++vx.i)
			v(vx.i, vy.i) = V(vx, vy, vx.i, vy.i);

	for(py.i = -(int)margin; py.i < (int)(p.Ny + margin); ++py.i)
		for(px.i = -(int)margin; px.i < (int)(p.Nx + margin); ++px.i)
			p(px.i, py.i) = P(px, py, px.i, py.i);
}
