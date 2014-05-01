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
	g_arrow(x0, y0, x1, y1, 0.04*sqrt(x1*x1 + y1*y1), .004);
}

int main()
{
	const double Re = 50;//Reynolds
	const double a = 1;//parabola parameter
	const unsigned int N = 8;//x-y divide num
	const double width = 1;//length from origin to bound
	const double arrow_len = 2*width/(N*a);

	//NS_System nssys(Re, a, N, width);
	NS_Staggered nssys(Re, N, N, [](double x, double y)->double{return 1;}, [](double x, double y)->double{return 1;});

	g_init("", 200, 200);
	g_device(G_DISP);
	g_def_scale(0, -width * 1.1, width*1.1, -width*1.1, width*1.1, 0, 0, 200, 200);
	g_sel_scale(0);

	g_text_color(G_BLACK);
	g_text_font(G_FONT_TIMES_24);

//		for(nssys.j = 0; nssys.j < N; ++nssys.j)
//		{
//			for(nssys.i = 0; nssys.i < N; ++nssys.i)
//				std::cout << nssys.u.u[nssys.n] << "\t" << nssys.u.v[nssys.n] << std::endl;
//			std::cout << std::endl;
//		}

	while(1)
	{
		g_cls();

		g_line_color(G_BLACK);
		g_move(-width, width);
		g_plot(-width, -width);
		g_plot(width, -width);
		g_plot(width, width);

		for(nssys.j = 0; nssys.j < N + 1; ++nssys.j)
			for(nssys.i = 0; nssys.i < N + 1; ++nssys.i)
				g_arrow2(nssys.px, nssys.py, nssys.u.u[nssys.un], nssys.u.v[nssys.un]);
		
		g_sleep(-1);

		std::cout << nssys.ux << std::endl;
	}

	g_term();
	return 0;
}
