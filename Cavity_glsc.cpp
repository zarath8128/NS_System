#include "NS.h"
#include <ZNAC/LA/LEQSolver.h>
#include <iostream>
#include <glsc.h>
#include <ctime>
#include <cmath>
#include <cstdint>

using namespace NS;
using namespace ZNAC::LA;

void g_arrow2(double x0, double y0, double x1, double y1, double len)
{
	g_arrow(x0, y0, x1, y1, 0.6*len*sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0)), .01);
}



void option_analays(int argc, char *argv[]);

/*typedef union D
{
	double d;
	uint64_t n;
}D;*/

void printval(int Nx, int Ny, NS_Vector &v, const char *description, bool Negative)
{
	D d;
	std::cout << description << std::endl;
	std::cout << "u:\n";
	for(int j = 0; j < Ny; ++j)
	{
		for(int i = 0; i < Nx; ++i)
		{
			if(Negative)
				d.d = -v.u[j*Nx + Nx - 1 - i];
			else
				d.d = v.u[j*Nx + i];
			std::cout << std::hex << (d.d == 0?0:d.n) << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\nv:\n";
	for(int j = 0; j < Ny; ++j)
	{
		for(int i = 0; i < Nx; ++i)
		{
			if(Negative)
				d.d = v.v[j*Nx + Nx - 1 - i];
			else
				d.d = v.v[j*Nx + i];
			std::cout << std::hex << (d.d == 0?0:d.n) << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";
}

int main(int argc, char *argv[])
{
//	NS_Option opt(argc, argv);
	double Re = 400;//Reynolds
	double a = 1;//parabola parameter
	unsigned int N = 64;//x-y divide num
	float width = 1;//length from origin to bound
	double arrow_len = 2*width/(16*a);
	bool Negative = true;
	int pn_time = 1;
	bool All = false;
	bool Advection = !All && true;
	bool Diffusion = !All && true;
	bool Force = !All && false;
	bool Pressure = !All && true;
	bool Auto = true;
	double dt = (2*width/N)*(2*width/N)*0.25*0.8*Re;
	double t = 0;
	float wnd_width = 300;
	const float bound_x[] = {-width, -width, width, width};
	const float bound_y[] = {width, -width, -width, width};
	char dammy_name[] = "";
	char info[256];
	int cut = 100;
	int count = -1;
	int pmod = N / 32;
	int vmod = N / 16;
	double draw_p[N/pmod][N/pmod];

	if(argc == 2)
	{
		sscanf(argv[1], "%lf", &Re);
	}

	a *= Negative?-1:1;

	//std::cout << "dt = " << dt << std::endl;

	NS_System nssys(Re, a, N, width);
	NS_Vector F(nssys.Nx, nssys.Ny);
	NS_Diffusion diffusion(nssys.Nx, nssys.Ny, nssys.Re, nssys.dx);
	NS_Diffusion_S diffusion_s(nssys.Nx, nssys.Ny, nssys.dx);
	NS_Advection advection(nssys.Nx, nssys.Ny, nssys.dx);
	NS_Vector vtmp1(nssys.Nx, nssys.Ny);
	Vector<double> stmp(nssys.Nx*nssys.Ny), phi(nssys.Nx*nssys.Ny);


	for(unsigned int j = 0; j < N; ++j)
		for(unsigned int i = 0; i < N; ++i)
			F.u[j*nssys.Nx + i] = F.v[j*nssys.Nx + i] = 0;

	for(unsigned int j = 0; j < N; ++j)
		for(unsigned int i = 0; i < N; ++i)
		{
			vtmp1.u[j*nssys.Nx + i] = 0;
			vtmp1.v[j*nssys.Nx + i] = 0;
			phi[j*nssys.Nx + i] = 0;
			stmp[j*nssys.Nx + i] = 0;
		}
		

	g_init(dammy_name, wnd_width, wnd_width);
	g_device(G_DISP);
	g_def_scale(0, -width * 1.1, width*1.1, -width*1.1, width*1.1, 0, 0, wnd_width, wnd_width);
	g_sel_scale(0);

	g_text_color(G_BLACK);
	g_text_font(G_FONT_TIMES_24);

	g_capture_set("");

	while(1)
	{
		std::cerr << count << "\n";
		/*for(unsigned int i = 0; i < N; ++i)
			if(nssys.x(i) != -nssys.x(nssys.Nx - 1 - i))
			{
				D d1, d2;
				d1.d = nssys.x(i);
				d2.d = -nssys.x(nssys.Nx - 1 - i);
				std::cout << "x" << i << ":" << std::endl;
				std::cout << std::hex << d1.n << std::endl;
				std::cout << std::hex << d2.n << std::endl;
				break;
			}

		for(unsigned int i = 0; i < N; ++i)
			if(nssys.y(i) != -nssys.y(nssys.Ny - 1 - i))
			{
				D d1, d2;
				d1.d = nssys.y(i);
				d2.d = -nssys.y(nssys.Ny - 1 - i);
				std::cout << "y" << i << ":" << std::endl;
				std::cout << std::hex << d1.n << std::endl;
				std::cout << std::hex << d2.n << std::endl;
				break;
			}
*/

		if(count ++ % cut == 0)
		{
		if(Auto)
			g_sleep(0.016);
		else
			g_sleep(-1);
			g_cls();

			for(unsigned int j = 0; j < N; ++j)
				for(unsigned int i = 0; i < N; ++i)
					if(i % pmod == 0 && j % pmod == 0)
						draw_p[i/pmod][j/pmod] = nssys.p[j*N + i];

			g_line_color(G_RED);
			for(double lev = 0; lev < 40; lev += 0.02)
				g_contln(-width , width , -width , width , (G_REAL *)draw_p, N/pmod, N/pmod, lev);

			g_line_color(G_BLACK);
			for(unsigned int j = 1; j < N; ++j)
				for(unsigned int i = 1; i < N - 1; ++i)
					if((i + 1)% vmod == 0 && (j + 1) % vmod == 0)
						g_arrow2(nssys.x(i), nssys.y(j), nssys.u.u[j*N + i], nssys.u.v[j*N + i], arrow_len);
	
			g_line_color(G_BLUE);
			g_move(-width, width);
			g_plot(-width, -width);
			g_plot(width, -width);
			g_plot(width, width);

			sprintf(info, "t:%f  Re:%f", t, Re);
			g_text(4, 10, info);
			g_capture();
		}

		if(pn_time == count)
			printval(N, N, nssys.u , "p-Gradient", Negative);
		
		if(All)
		{
			NS_Vector diff(nssys.Nx, nssys.Ny), adv(nssys.Nx, nssys.Ny);

			advection(nssys.u, adv);
			diffusion(nssys.u, diff);
		}

		if(Advection)
		{
			advection(nssys.u, vtmp1);
			if(pn_time == count)
				printval(N, N, vtmp1, "Advection", Negative);
			for(unsigned int j = 1; j < N - 1; ++j)
				for(unsigned int i = 1; i < N - 1; ++i)
				{	
					nssys.u.u[j*nssys.Nx + i] -= dt*vtmp1.u[j*nssys.Nx + i];
					nssys.u.v[j*nssys.Nx + i] -= dt*vtmp1.v[j*nssys.Nx + i];
				}
		}

		if(Diffusion)
		{
			diffusion(nssys.u, vtmp1);
			if(pn_time == count)
				printval(N, N, vtmp1, "Diffusion", Negative);
			for(unsigned int j = 1; j < N - 1; ++j)
				for(unsigned int i = 1; i < N - 1; ++i)
				{
					nssys.u.u[j*nssys.Nx + i] += dt*vtmp1.u[j*nssys.Nx + i];
					nssys.u.v[j*nssys.Nx + i] += dt*vtmp1.v[j*nssys.Nx + i];
				}
		}

		if(Force)
		{
			for(unsigned int j = 1; j < N - 1; ++j)
				for(unsigned int i = 1; i < N - 1; ++i)
				{
					nssys.u.u[j*nssys.Nx + i] += dt*F.u[j*nssys.Nx + i];
					nssys.u.v[j*nssys.Nx + i] += dt*F.v[j*nssys.Nx + i];
				}
		}

		if(Pressure)
		{
			NS_Gradient(nssys.p, vtmp1, nssys.dx);
			if(pn_time == count)
				printval(N, N, vtmp1, "p-Gradient", Negative);
			for(unsigned int j = 1; j < N - 1; ++j)
				for(unsigned int i = 1; i < N - 1; ++i)
				{
					nssys.u.u[j*nssys.Nx + i] -= dt*vtmp1.u[j*nssys.Nx + i];
					nssys.u.v[j*nssys.Nx + i] -= dt*vtmp1.v[j*nssys.Nx + i];
				}
			NS_Divergence(nssys.u, stmp, nssys.dx);
			CG<double>(N*N, 1e-15, supNorm<double>())(diffusion_s, phi, stmp);

			NS_Gradient(phi, vtmp1, nssys.dx);
			if(pn_time == count)
				printval(N, N, vtmp1, "phi-Gradient", Negative);
			for(unsigned int j = 1; j < N - 1; ++j)
				for(unsigned int i = 1; i < N - 1; ++i)
				{
					nssys.u.u[j*nssys.Nx + i] += vtmp1.u[j*nssys.Nx + i];
					nssys.u.v[j*nssys.Nx + i] += vtmp1.v[j*nssys.Nx + i];
				}
			for(unsigned int j = 1; j < N - 1; ++j)
				for(unsigned int i = 1; i < N - 1; ++i)
					nssys.p[j*nssys.Nx + i] -= 1/dt * phi[j*nssys.Nx + i];

			for(unsigned int i = 0; i < N; ++i)
				nssys.p[i] = nssys.p[nssys.Nx + i], nssys.p[(nssys.Ny - 1)*nssys.Nx + i] = nssys.p[(nssys.Ny - 2)*nssys.Nx + i];
			for(unsigned int i = 0; i < N; ++i)
				nssys.p[i*nssys.Nx] = nssys.p[i*nssys.Nx + 1], nssys.p[i*nssys.Nx + nssys.Nx - 1] = nssys.p[i*nssys.Nx + nssys.Nx - 2];

			double pmin = 0, pmax = 0;
			for(unsigned int i = 0; i < N*N; ++i)
				pmin = (pmin < nssys.p[i] ? pmin : nssys.p[i]), pmax = (pmax > nssys.p[i] ? pmax : nssys.p[i]);
			for(unsigned int i = 0; i < N*N; ++i)
				nssys.p[i] -=  pmin;
		}

		t += dt;

		if(pn_time == count)
			break;

	}

	g_term();
	return 0;
}
