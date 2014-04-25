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
	g_arrow(x0, y0, x1, y1, 0.01*sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0)), .001);
}

int main()
{
	const double Re = 50;//Reynolds
	const double a = 1;//parabola parameter
	const unsigned int N = 16;//x-y divide num
	const double width = 1;//length from origin to bound
	const double arrow_len = 2*width/(N*a);

	NS_System nssys(Re, a, N, width);

	g_init("", 200, 200);
	g_device(G_DISP);
	g_def_scale(0, -width * 1.1, width*1.1, -width*1.1, width*1.1, 0, 0, 200, 200);
	g_sel_scale(0);

	g_text_color(G_BLACK);
	g_text_font(G_FONT_TIMES_24);

	while(1)
	{
		g_cls();

		for(unsigned int j = 0; j < N; ++j)
			for(unsigned int i = 0; i < N; ++i)
				g_arrow2(nssys.x(i), nssys.y(j), nssys.u.u[j*N + i], nssys.u.v[j*N + i]);
		
		g_sleep(0.016);
	}

	g_term();
	return 0;
}
