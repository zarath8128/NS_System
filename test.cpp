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
	g_arrow(x0, y0, x1, y1, 0.08*sqrt(x1*x1 + y1*y1), .01);
}

double top(double x){return 1;}
void boundary(NS_Staggered &ns);

int main()
{
	const double Re = 50;//Reynolds
	const double a = 1;//parabola parameter
	const unsigned int N = 128;//x-y divide num
	const double width = 1;//length from origin to bound
	const double dx = 2*width/N;
	const double dy = 2*width/N;
	char dammy_name[] = "";

	const bool WriteOut = false;
	const bool Auto = true;
	const bool Skip = false;
	const bool Diffusion = true;
	const bool Advection = true;
	const bool Pressure = true;

	//auto top = [a](double x)->double{return a*(1-x)*(1+x);};

	NS_Staggered nssys(
			Re, N, N 
			//[](double x, double y, int i, int j)->double{return j == N - 2?5:0;});
			//[](double x, double y, int i, int j)->double{return 0;}
	);

	boundary(nssys);

	NS_Grid utmp(N + 1, N), vtmp(N, N + 1), phi(N, N), stmp(N, N);
	for(unsigned int i = 0; i < phi.N(); ++i)
		phi[i] = stmp[i] = 0;

	NS::Diffusion diffusion_u(1/(Re*dx*dx), 1/(Re*dy*dy), N + 1, N);
	NS::Diffusion diffusion_v(1/(Re*dx*dx), 1/(Re*dy*dy), N, N + 1);
	NS::Diffusion diffusion_p(-1/(dx*dx), -1/(dy*dy), N, N);
	NS::Advection advection(dx, dy);
	NS::Gradient  gradient(dx, dy);
	NS::Divergence divergence(dx, dy);

	const double dt = 0.8*0.5*Re/(1/(dx*dx) + 1/(dy*dy));

	g_init(dammy_name, 200, 200);
	g_device(G_DISP);
	g_def_scale(0, -width * 1.1, width*1.1, -width*1.1, width*1.1, 0, 0, 200, 200);
	g_sel_scale(0);

	g_text_color(G_BLACK);
	g_text_font(G_FONT_TIMES_24);

	while(1)
	{
		g_cls();

		g_line_color(G_BLACK);
		g_move(-width, width);
		g_plot(-width, -width);
		g_plot(width, -width);
		g_plot(width, width);

		g_arrow2(-width*1.05, width*1.01, 1, 0);
		g_arrow2(-width*1.05, width*1.01, 0, 1);

		{
			g_line_color(G_RED);
			double dl = 0.1;
			for(double lev = -10; lev < 10; lev += dl)
				g_contln(-width, width, -width, width, (G_REAL *)(double*)nssys.p, N + 2, N + 2, lev);
		}
		/*{
			int &i = nssys.px.i;
			int &j = nssys.py.i;
			for(j = 0; j < nssys.p.Ny; ++j)
				for(i = 0; i < nssys.p.Nx; ++i)
					g_arrow2(nssys.px, nssys.py, (nssys.u(i, j) + nssys.u(i + 1, j))*0.5, (nssys.v(i, j) + nssys.v(i, j + 1))*0.5);
		}*/

		if(Auto)
			if(!Skip)
				g_sleep(0.16);
			else;
		else
			g_sleep(-1);

		if(Diffusion)
		{
			diffusion_u(nssys.u, utmp);
			diffusion_v(nssys.v, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					nssys.u(xi, yi) += dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.v(xi, yi) += dt*vtmp(xi, yi);
			boundary(nssys);
		}

		if(Advection)
		{
			advection(nssys.u, nssys.v, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					nssys.u(xi, yi) += dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.v(xi, yi) += dt*vtmp(xi, yi);
			boundary(nssys);
		}

		if(Pressure)
		{
			gradient(nssys.p, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					nssys.u(xi, yi) -= dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.v(xi, yi) -= dt*vtmp(xi, yi);
			//boundary(nssys);


			divergence(nssys.u, nssys.v, stmp);

			/*for(int i = -1; i < (int)N + 1; ++i)
			{
				for(int j = -1; j < (int)N + 1; ++j)
					std::cout << stmp(j, i) << "\t";
				std::cout << std::endl;
			}*/
			CG<double>(N*N, 1e-10, supNorm<double>())(diffusion_p, phi, stmp);

			gradient(phi, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					nssys.u(xi, yi) += dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.v(xi, yi) += dt*vtmp(xi, yi);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.p(xi, yi) -= phi(xi, yi)/dt;
			boundary(nssys);
		}

	}

	g_term();
	return 0;
}

void boundary(NS_Staggered &ns)
{
	{
		int &i = ns.ux.i;
		int Nx = ns.u.Nx;
	       	int margin = ns.u.margin;
		for(i = -margin; i < Nx + margin ; ++i)
		{
			ns.u(i, ns.u.Ny) = 2*top(ns.ux) - ns.u(i, ns.u.Ny - 1);
			ns.u(i, -1) = - ns.u(i, 0);
		}
	}

	{
		int &j = ns.uy.i;
		int Ny = ns.u.Ny;
	       	int margin = ns.u.margin;
		for(j = -margin; j < Ny + margin ; ++j)
		{
			ns.u(ns.u.Nx - 1, j) = 0;
			ns.u(0, j) = 0;
		}
	}

	{
		int &i = ns.vx.i;
		int Nx = ns.v.Nx;
	       	int margin = ns.v.margin;
		for(i = -margin; i < Nx + margin ; ++i)
		{
			ns.v(i, ns.v.Ny - 1) = 0;
			ns.v(i, 0) = 0;
		}
	}

	{
		int &j = ns.vx.i;
		int Ny = ns.v.Ny;
	       	int margin = ns.v.margin;
		for(j = -margin; j < Ny + margin ; ++j)
		{
			ns.v(ns.v.Nx, j) = - ns.v(ns.v.Nx - 1, j);
			ns.v(-1, j) = - ns.v(0, j);
		}
	}

	{
		int &i = ns.px.i;
		int Nx = ns.p.Nx;
	       	int margin = ns.p.margin;
		for(i = -margin; i < Nx + margin ; ++i)
		{
			ns.p(i, ns.p.Ny) = ns.p(i, ns.p.Ny - 1);
			ns.p(i, -1) = ns.p(i, 0);
		}
	}

	{
		int &j = ns.py.i;
		int Ny = ns.p.Ny;
	       	int margin = ns.p.margin;
		for(j = -margin; j < Ny + margin ; ++j)
		{
			ns.p(ns.p.Nx, j) = ns.p(ns.p.Nx - 1, j);
			ns.p(-1, j) = ns.p(0, j);
		}
	}
}
