#include "NS.h"
#include <ZNAC/LA/LEQSolver.h>
#include <iostream>
#include <glsc.h>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <Rainbow.c>

#define SWITCHTEXT(b) (b?"On":"Off")

using namespace NS;
using namespace ZNAC::LA;

void g_arrow2(double x0, double y0, double x1, double y1)
{
	g_arrow(x0, y0, x1, y1, 0.08*sqrt(x1*x1 + y1*y1), .01);
}

constexpr double a = 0.5, b = 0;

double top(double x){return a*(1 + x)*(1 - x) + b;}
//double top(double x){return a;}
void boundary(NS_Staggered &ns);

int main()
{
	const double Re = 400;//Reynolds
	const unsigned int N = 64;//x-y divide num
	const double width = 1;//length from origin to bound
	const double dx = 2*width/N;
	const double dy = 2*width/N;
	const double cg_err = 1e-6;
	double t = 0;
	const double dt = 0.8*0.5*Re/(1/(dx*dx) + 1/(dy*dy));
	char dammy_name[] = "";
	const unsigned int ContNum = 200;
	const unsigned int Cut = 10;
	unsigned int counter = 0;
	double div_max = 0, div_min = 0;

	const bool WriteOut = false;
	const bool Auto = true;
	const bool Skip = false;
	const bool Diffusion = true;
	const bool Advection = true;
	const bool Pressure = true;
	const bool All = false;
	const bool P_Grid = false;
	const bool U_Grid = false;
	const bool V_Grid = false;

	NS_Staggered nssys(
			Re, N, N 
			//[](double x, double y, int i, int j)->double{return j == N - 2?5:0;});
			//[](double x, double y, int i, int j)->double{return 0;}
	);

	boundary(nssys);

	NS_Grid utmp(N + 1, N), vtmp(N, N + 1), utmpa(N + 1, N), vtmpa(N, N + 1), phi(N, N), stmp(N, N);

	NS::Diffusion diffusion_u(1/(Re*dx*dx), 1/(Re*dy*dy), N + 1, N);
	NS::Diffusion diffusion_v(1/(Re*dx*dx), 1/(Re*dy*dy), N, N + 1);
	NS::Diffusion diffusion_p(-1/(dx*dx), -1/(dy*dy), N, N);
	NS::Advection advection(dx, dy);
	NS::Gradient  gradient(dx, dy);
	NS::Divergence divergence(dx, dy);
	CG<double> cg(N*N, cg_err, supNorm<double>());


	g_init(dammy_name, 300, 200);
	g_device(G_DISP);
	g_def_scale(0, -width * 1.1, width*1.1, -width*1.1, width*1.1, 0, 0, 200, 200);
	g_sel_scale(0);

	g_text_color(G_BLACK);
	g_text_font(G_FONT_TIMES_24);

	while(1)
	{
		if(counter ++ % Cut == 0){
		g_cls();

		g_line_color(G_BLACK);
		g_move(-width, width);
		g_plot(-width, -width); g_plot(width, -width);
		g_plot(width, width);

		g_arrow2(-width*1.05, width*1.01, 1, 0);
		g_arrow2(-width*1.05, width*1.01, 0, 1);

		{
			double P[N][N];
			double min = DBL_MAX, max = DBL_MIN; 
			for(int i = 0; i < N; ++i)
				for(int j = 0; j < N; ++j)
				{
					min = min > nssys.p(i, j) ? nssys.p(i, j) : min;
					max = max < nssys.p(i, j) ? nssys.p(i, j) : max;
				}

			max -= min;

			for(int i = 0; i < N; ++i)
				for(int j = 0; j < N; ++j)
					P[i][j] = (nssys.p(i, j) -= min)/max;

			const double d = 1./(ContNum - 1);
			double r = 1, g = 0, b = 0;
			for(int i = 0; i < ContNum; ++i)
			{
				Rainbow(1 - 1.*d*i, 1, 1, &r, &g, &b);
				g_line_color(g_rgb_color(r, g, b));
				g_contln(-width + 0.5*dx, width - 0.5*dx, -width + 0.5*dx, width - 0.5*dx, (G_REAL *)(double*)P, N, N, d*i);
			}
		}
		if(P_Grid)
		{
			g_line_color(G_GREEN);
			int &i = nssys.ux.i;
			int &j = nssys.vy.i;
			for(j = 0; j < (int)(nssys.p.Ny + 1); ++j)
				g_move(-width, nssys.vy), g_plot(width, nssys.vy);
			for(i = 0; i < (int)(nssys.p.Nx + 1); ++i)
				g_move(nssys.ux, -width), g_plot(nssys.ux, width);
					
		}

		if(U_Grid)
		{
			g_line_color(G_BLUE);
			int &i = nssys.px.i;
			int &j = nssys.vy.i;
			for(j = 0; j < (int)(nssys.p.Ny + 1); ++j)
				g_move(-width - 0.5*dx, nssys.vy), g_plot(width + 0.5*dx, nssys.vy);
			for(i = -1; i < (int)(nssys.p.Nx + 1); ++i)
				g_move(nssys.px, -width), g_plot(nssys.px, width);
					
		}

		if(V_Grid)
		{
			g_line_color(G_MAGENTA);
			int &i = nssys.ux.i;
			int &j = nssys.py.i;
			for(j = -1; j < (int)(nssys.p.Ny + 1); ++j)
				g_move(-width, nssys.py), g_plot(width, nssys.py);
			for(i = 0; i < (int)(nssys.p.Nx + 1); ++i)
				g_move(nssys.ux, -width - 0.5*dy), g_plot(nssys.ux, width + 0.5*dy);
					
		}

		{
			g_line_color(G_BLACK);
			int &i = nssys.px.i;
			int &j = nssys.py.i;
			for(j = 0; j < nssys.p.Ny; j += 1)
				for(i = 0; i < nssys.p.Nx; i += 1)
					g_arrow2(nssys.px, nssys.py, (nssys.u(i, j) + nssys.u(i + 1, j))*0.5, (nssys.v(i, j) + nssys.v(i, j + 1))*0.5);
		}

		/*---text block--*/
		char text[1024];
		sprintf(text , "t = %5.2f, dt = %5.2f", t, dt);
		g_text(200.0, 10, text);
		sprintf(text , "Re = %5.2f", Re);
		g_text(200.0, 20, text);
		sprintf(text , "a = %5.2f, b = %5.2f", a, b);
		g_text(200.0, 30, text);
		sprintf(text , "dx = %5.2e, dy = %5.2e", dx, dy);
		g_text(200.0, 40, text);
		sprintf(text , "N = %d", N);
		g_text(200.0, 50, text);
		sprintf(text, "Div. Max = %5.2e", div_max);
		g_text(200.0, 60, text);
		sprintf(text, "Div. Min = %5.2e", div_min);
		g_text(200.0, 70, text);
		sprintf(text , "Diffusion:%s", SWITCHTEXT(Diffusion));
		g_text(200.0, 80, text);
		sprintf(text , "Advection:%s", SWITCHTEXT(Advection));
		g_text(200.0, 90, text);
		sprintf(text , "Pressure :%s", SWITCHTEXT(Pressure));
		g_text(200.0, 100, text);
		sprintf(text , "Auto     :%s", SWITCHTEXT(Auto));
		g_text(200.0, 110, text);
		sprintf(text , "Skip     :%s", SWITCHTEXT(Skip));
		g_text(200.0, 120, text);
		sprintf(text , "Cut      :%d", Cut);
		g_text(200.0, 130, text);
		sprintf(text , "CG error limit:%e", cg_err);
		g_text(200.0, 140, text);
		sprintf(text , "CG loop limit:%d", N*N);
		g_text(200.0, 150, text);
		sprintf(text , "CG loos:%d", cg.n);
		g_text(200.0, 160, text);


		}

		if(Auto)
			if(!Skip)
				g_sleep(0.16);
			else;
		else
			g_sleep(-1);

		if(Diffusion && !All)
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

		if(Advection && !All)
		{
			advection(nssys.u, nssys.v, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					nssys.u(xi, yi) -= dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.v(xi, yi) -= dt*vtmp(xi, yi);
			boundary(nssys);
		}

		if(Pressure && !All)
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
			cg(diffusion_p, phi, stmp);

			gradient(phi, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					nssys.u(xi, yi) += utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.v(xi, yi) += vtmp(xi, yi);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.p(xi, yi) -= phi(xi, yi)/dt;
		
			divergence(nssys.u, nssys.v, stmp);

			div_max = stmp(0, 0), div_min = stmp(0, 0);
			for(int i = 0; i < N; ++i)
				for(int j = 0; j < N; ++j)
				{
					div_max = div_max < stmp(i, j) ? stmp(i, j) : div_max;
					div_min = div_min > stmp(i, j) ? stmp(i, j) : div_min;
				}
			//std::cout << "max : " << max << std::endl;
			//std::cout << "min : " << min << std::endl;

				boundary(nssys);
		}

		if(All)
		{
			diffusion_u(nssys.u, utmp);
			diffusion_v(nssys.v, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					utmpa(xi, yi) = nssys.u(xi, yi) + dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					vtmpa(xi, yi) = nssys.v(xi, yi) + dt*vtmp(xi, yi);
			
			advection(nssys.u, nssys.v, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					utmpa(xi, yi) -= dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					vtmpa(xi, yi) -= dt*vtmp(xi, yi);

			gradient(nssys.p, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					utmpa(xi, yi) -= dt*utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					vtmpa(xi, yi) -= dt*vtmp(xi, yi);

			divergence(utmpa, vtmpa, stmp);

			CG<double>(-1, 1e-10, supNorm<double>())(diffusion_p, phi, stmp);

			gradient(phi, utmp, vtmp);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N + 1; ++xi)
					nssys.u(xi, yi) = utmpa(xi, yi) + utmp(xi, yi);
			for(int yi = 0; yi < N + 1; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.v(xi, yi) = vtmpa(xi, yi) + vtmp(xi, yi);
			for(int yi = 0; yi < N; ++yi)
				for(int xi = 0; xi < N; ++xi)
					nssys.p(xi, yi) -= phi(xi, yi)/dt;
			boundary(nssys);
		}
		t += dt;
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
		for(i = 0; i < Nx; ++i)
		{
			ns.u(i, ns.u.Ny) = 2*top(ns.ux) - ns.u(i, ns.u.Ny - 1);
			ns.u(i, -1) = - ns.u(i, 0);
		}
	}

	{
		int &j = ns.uy.i;
		int Ny = ns.u.Ny;
	       	int margin = ns.u.margin;
		for(j = 0; j < Ny; ++j)
		{
			ns.u(ns.u.Nx - 1, j) = 0;
			ns.u(0, j) = 0;
		}
	}

	{
		int &i = ns.vx.i;
		int Nx = ns.v.Nx;
	       	int margin = ns.v.margin;
		for(i = 0; i < Nx; ++i)
		{
			ns.v(i, ns.v.Ny - 1) = 0;
			ns.v(i, 0) = 0;
		}
	}

	{
		int &j = ns.vy.i;
		int Ny = ns.v.Ny;
	       	int margin = ns.v.margin;
		for(j = 0; j < Ny; ++j)
		{
			ns.v(ns.v.Nx, j) = - ns.v(ns.v.Nx - 1, j);
			ns.v(-1, j) = - ns.v(0, j);
		}
	}

	{
		int &i = ns.px.i;
		int Nx = ns.p.Nx;
	       	int margin = ns.p.margin;
		for(i = 0; i < Nx ; ++i)
		{
			ns.p(i, ns.p.Ny) = ns.p(i, ns.p.Ny - 1) - 1/(ns.Re*ns.dy)*(2*ns.p(i, ns.p.Ny - 1) - 5*ns.p(i, ns.p.Ny - 2) + 4*ns.p(i, ns.p.Ny - 3) - ns.p(i, ns.p.Ny - 4));
			ns.p(i, -1) = ns.p(i, 0) + 1/(ns.Re*ns.dy)*(2*ns.p(i, 0) - 5*ns.p(i, 1) + 4*ns.p(i, 2) - ns.p(i, 3));
		}
	}

	{
		int &j = ns.py.i;
		int Ny = ns.p.Ny;
	       	int margin = ns.p.margin;
		for(j = 0; j < Ny; ++j)
		{
			ns.p(ns.p.Nx, j) = ns.p(ns.p.Nx - 1, j) - 1/(ns.Re*ns.dx)*(2*ns.p(ns.p.Nx - 1, j) - 5*ns.p(ns.p.Nx - 2, j) + 4*ns.p(ns.p.Nx - 3, j) - ns.p(ns.p.Nx - 4, j));
			ns.p(-1, j) = ns.p(0, j) + 1/(ns.Re*ns.dx)*(2*ns.p(0, j) - 5*ns.p(0, j) + 4*ns.p(0, j) - ns.p(0, j));
		}
	}
}
