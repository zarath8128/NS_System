#include "NS.h"
#include <ZNAC/LA/LEQSolver.h>
#include <iostream>
#include <glsc.h>
#include <ctime>
#include <cmath>

using namespace NS;
using namespace ZNAC::LA;

void g_arrow2(double x0, double y0, double x1, double y1)
{
	g_arrow(x0, y0, x1, y1, 0.1*sqrt(x1*x1 + y1*y1), .01);
}

int main()
{
	const double Re = 10;//Reynolds
	const double a = 1;//parabola parameter
	const unsigned int N = 16;//x-y divide num
	const double width = 1;//length from origin to bound
	const double arrow_len = 2*width/(N*a);
	char dammy_name[] = "";

	const bool WriteOut = false;
	const bool Auto = true;
	const bool Diffusion = true;
	const bool Advection = true;
	const bool Pressure = true;

	//auto top = [a](double x)->double{return a*(1-x)*(1+x);};
	auto top = [a](double x)->double{return 1;};

	NS_Staggered nssys(Re, N, N);
	NS_Vector vtmp(N + 2, N + 3);
	ZNAC::LA::Vector<double> phi((N + 2)*(N + 2));
	ZNAC::LA::Vector<double> stmp((N + 2)*(N + 2));
	for(unsigned int i = 0; i < phi.N(); ++i)
		phi[i] = stmp[i] = 0;

	int &i = nssys.i, &j = nssys.j;

	const double dt = 0.8*0.5*nssys.Re/(1/(nssys.dx*nssys.dx) + 1/(nssys.dy*nssys.dy));

	g_init(dammy_name, 200, 200);
	g_device(G_DISP);
	g_def_scale(0, -width * 1.1, width*1.1, -width*1.1, width*1.1, 0, 0, 200, 200);
	g_sel_scale(0);

	g_text_color(G_BLACK);
	g_text_font(G_FONT_TIMES_24);

	while(1)
	{
		//std::cout << "dt:" << dt << std::endl;
		//std::cout << "redx2_inv:" << nssys.redx2_inv << std::endl;
		//std::cout << "dx_inv:" << nssys.dx_inv << std::endl;
		
		if(WriteOut)
		{
			std::cout << "======================" << std::endl;
			for(j = 0; j < nssys.Ny + 1; ++j)
			{
				for(i = 0; i < nssys.Nx + 1; ++i)
					std::cout << nssys.un << ":" << nssys.u.u[nssys.un] << std::endl;
				std::cout << std::endl;
			}
		}

		g_cls();

		g_line_color(G_BLACK);
		g_move(-width, width);
		g_plot(-width, -width);
		g_plot(width, -width);
		g_plot(width, width);

		g_line_color(G_RED);
		for(i = 0; i < nssys.Nx; ++i)
			g_arrow2(nssys.vx, 1, top(nssys.vx), 0);

		g_line_color(G_BLACK);
		for(j = 0; j < nssys.Ny; ++j)
			for(i = 0; i < nssys.Nx; ++i)
				g_arrow2(nssys.px, nssys.py, (nssys.u.u[nssys.un] + nssys.u.u[nssys.un(i + 1, j)])*.5, (nssys.u.v[nssys.vn] + nssys.u.v[nssys.vn(i, j + 1)])*.5);
		
		if(Auto)
			g_sleep(0.016);
		else
			g_sleep(-1);

		if(Diffusion)
		{
			nssys.Diffusion_u(vtmp);
			j = nssys.Ny - 1;
			for(i = 1; i < nssys.Nx; ++i)
				vtmp.u[nssys.un] = nssys.redx2_inv*(2*top(nssys.ux) + nssys.u.u[nssys.un(i, j - 1)] - 3*nssys.u.u[nssys.un])
					+ nssys.redy2_inv*(nssys.u.u[nssys.un(i + 1, j)] + nssys.u.u[nssys.un(i - 1, j)] - 2*nssys.u.u[nssys.un]);
			j = nssys.Ny;
			for(i = 0; i < nssys.Nx; ++i)
				vtmp.v[nssys.vn] = 0;

			j = 0;
			for(i = 1; i < nssys.Nx; ++i)
				vtmp.u[nssys.un] = nssys.redx2_inv*(nssys.u.u[nssys.un(i, j + 1)] - 3*nssys.u.u[nssys.un])
					+ nssys.redy2_inv*(nssys.u.u[nssys.un(i + 1, j)] + nssys.u.u[nssys.un(i - 1, j)] - 2*nssys.u.u[nssys.un]);
			for(i = 0; i < nssys.Nx; ++i)
				vtmp.v[nssys.vn] = 0;

			i = 0;
			for(j = 1; j < nssys.Ny; ++j)
			{
				vtmp.u[nssys.un] = 0;
				vtmp.v[nssys.vn] = nssys.redx2_inv*(nssys.u.v[nssys.vn(i + 1, j)] - 3*nssys.u.v[nssys.vn])
					+ nssys.redy2_inv*(nssys.u.v[nssys.vn(i, j + 1)] + nssys.u.v[nssys.vn(i - 1, j)] - 2*nssys.u.v[nssys.vn]);
			}
			i = nssys.Nx;
			for(j = 1; j < nssys.Ny; ++j)
				vtmp.u[nssys.un] = 0;
			i = nssys.Nx - 1;
			for(j = 1; j < nssys.Ny; ++j)
				vtmp.v[nssys.vn] = nssys.redx2_inv*(nssys.u.v[nssys.vn(i + 1, j)] - 3*nssys.u.v[nssys.vn])
					+ nssys.redy2_inv*(nssys.u.v[nssys.vn(i, j + 1)] + nssys.u.v[nssys.vn(i, j - 1)] - 2*nssys.u.v[nssys.vn]);

			for(j = 0; j < nssys.Ny; ++j)
				for(i = 0; i < nssys.Nx + 1; ++i)
					nssys.u.u[nssys.un] += dt*vtmp.u[nssys.un];
			for(j = 0; j < nssys.Ny + 1; ++j)
				for(i = 0; i < nssys.Nx; ++i)
					nssys.u.v[nssys.vn] += dt*vtmp.v[nssys.vn];
		}

		if(Advection)
		{
			nssys.Advection(vtmp);
			i = 0;
			for(j = 0; j < nssys.Ny; ++j)
				vtmp.u[nssys.un] = 0;
			i = nssys.Nx;
			for(j = 0; j < nssys.Ny; ++j)
				vtmp.u[nssys.un] = 0;
			j = nssys.Ny - 1;
			for(i = 1; i < nssys.Nx; ++i)
				vtmp.u[nssys.un] = nssys.u.v[nssys.vn]*nssys.dy_inv*(2*top(nssys.ux) -(nssys.u.u[nssys.un(i, j - 1)] + nssys.u.u[nssys.un]))
					+ nssys.u.u[nssys.un]*nssys.dx_inv*(nssys.u.u[nssys.un(i + 1, j)] - nssys.u.u[nssys.un(i - 1, j)]);
			j = nssys.Ny;
			for(i = 0; i < nssys.Nx; ++i)
				vtmp.v[nssys.vn] = 0;

			j = 0;
			for(i = 1; i < nssys.Nx; ++i)
				vtmp.u[nssys.un] = nssys.u.v[nssys.vn]*nssys.dy_inv*( -(nssys.u.u[nssys.un(i, j - 1)] + nssys.u.u[nssys.un]))
					+ nssys.u.u[nssys.un]*nssys.dx_inv*(nssys.u.u[nssys.un(i + 1, j)] - nssys.u.u[nssys.un(i - 1, j)]);
			for(i = 0; i < nssys.Nx; ++i)
				vtmp.v[nssys.vn] = 0;

			i = 0;
			for(j = 1; j < nssys.Ny; ++j)
				vtmp.v[nssys.vn] = nssys.u.v[nssys.vn]*nssys.dy_inv*( nssys.u.v[nssys.vn(i, j + 1)] - nssys.u.v[nssys.vn(i, j - 1)])
					+ nssys.u.u[nssys.un]*nssys.dx_inv*(-(nssys.u.v[nssys.vn(i + 1, j)] + nssys.u.v[nssys.vn]));

			i = nssys.Nx;
			for(j = 0; j < nssys.Ny; ++j)
				vtmp.u[nssys.un] = 0;
			i = nssys.Nx - 1;
			for(j = 1; j < nssys.Ny; ++j)
				vtmp.v[nssys.vn] = nssys.u.v[nssys.vn]*nssys.dy_inv*( nssys.u.v[nssys.vn(i, j + 1)] - nssys.u.v[nssys.vn(i, j - 1)])
					+ nssys.u.u[nssys.un]*nssys.dx_inv*(-(nssys.u.v[nssys.vn(i - 1, j)] + nssys.u.v[nssys.vn]));

			for(j = 0; j < nssys.Ny; ++j)
				for(i = 0; i < nssys.Nx + 1; ++i)
					nssys.u.u[nssys.un] -= dt*vtmp.u[nssys.un];
			for(j = 0; j < nssys.Ny + 1; ++j)
				for(i = 0; i < nssys.Nx; ++i)
					nssys.u.v[nssys.vn] -= dt*vtmp.v[nssys.vn];

		}

		if(Pressure)
		{
			nssys.Gradient(nssys.p, vtmp);
			for(j = 0; j < nssys.Ny; ++j)
				for(i = 0; i < nssys.Nx + 1; ++i)
					nssys.u.u[nssys.un] -= dt*vtmp.u[nssys.un];
			for(j = 0; j < nssys.Ny + 1; ++j)
				for(i = 0; i < nssys.Nx; ++i)
					nssys.u.v[nssys.vn] -= dt*vtmp.v[nssys.vn];
			nssys.Divergence(stmp);
			ZNAC::LA::CG<double>(N*N*N*N, 1e-10, supNorm<double>())(nssys.diffusion_p, phi, stmp);
			
			nssys.Gradient(phi, vtmp);

			for(j = 0; j < nssys.Ny; ++j)
				for(i = 0; i < nssys.Nx + 1; ++i)
					nssys.u.u[nssys.un] += vtmp.u[nssys.un];
			for(j = 0; j < nssys.Ny + 1; ++j)
				for(i = 0; i < nssys.Nx; ++i)
					nssys.u.v[nssys.vn] += vtmp.v[nssys.vn];

			for(unsigned int i = 0; i < phi.N(); ++i)
				nssys.p[i] -= phi[i]/dt;
			for(int i = 1; i < N + 1; ++i)
			{
				nssys.p[i] = nssys.p[N + 2 + i];
				nssys.p[(N + 2)*(N + 1) + i] = nssys.p[(N + 2)*N + i];
				nssys.p[i*(N + 2)] = nssys.p[i*(N + 2) + 1];
				nssys.p[i*(N + 2) + N + 1] = nssys.p[i*(N + 2) + N];
			}

		}
	}

	g_term();
	return 0;
}
