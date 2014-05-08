#ifndef ZARATH_NS_NS_H
#define ZARATH_NS_NS_H

#include <ZNAC/LA/Vector.h>
#include <ZNAC/LA/Matrix.h>
#include <cassert>

namespace NS
{
	class NS_Vector
		:public ZNAC::LA::IVector<double>
	{
	private:
		double *buf;

	public:
		const unsigned int Nx, Ny;
		double * const u, *const v;

		NS_Vector(unsigned int Nx, unsigned int Ny);
		~NS_Vector();

		double &operator[](unsigned int i);
		const double &operator[](unsigned int i) const;

		unsigned int N()const;

	};

	class NS_Matrix
		:public ZNAC::LA::IMatrix<double>
	{
	public:
		const unsigned int Nx, Ny;

		NS_Matrix(unsigned int Nx, unsigned int Ny);
	};

	class NS_Diffusion
		:public NS_Matrix
	{
	public:
		const double Redx2_inv;

		NS_Diffusion(unsigned int Nx, unsigned int Ny, double Re, double dx);
		double &operator()(unsigned int r, unsigned int c);
		const double &operator()(unsigned int r, unsigned int c)const;
		void operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const;
	};

	class NS_Diffusion_S
		:public NS_Matrix
	{
	public:
		const double dx2_inv;

		NS_Diffusion_S(unsigned int Nx, unsigned int Ny, double dx);
		double &operator()(unsigned int r, unsigned int c);
		const double &operator()(unsigned int r, unsigned int c)const;
		void operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const;
	};

	class NS_Advection
		:public NS_Matrix
	{
	public:
		const double dx_inv;

		NS_Advection(unsigned int Nx, unsigned int Ny, double dx);
		double &operator()(unsigned int r, unsigned int c);
		const double &operator()(unsigned int r, unsigned int c)const;
		void operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const;
	};

/*	template<class T>
	void NS_Gradient(const ZNAC::LA::IVector<T> &p, NS_Vector &cod, double dx);
	template<class T>
	void NS_Divergence(const NS_Vector &dom, ZNAC::LA::IVector<T> &cod, double dx);
	template<class T>
	void NS_Laplace(const ZNAC::LA::IVector<T> &dom, ZNAC::LA::IVector<T> &cod, unsigned int Nx, unsigned int Ny, double dx);
*/
	template<class T>
	void NS_Gradient(const ZNAC::LA::IVector<T> &p, NS_Vector &cod, double dx)
	{
		const double dx_inv = 1/(2*dx);
		for(unsigned int ny = 1; ny < cod.Ny - 1; ++ny)
			for(unsigned int nx = 1; nx < cod.Nx - 1; ++nx)
			{
				cod.u[ny*cod.Nx + nx] = dx_inv*(p[ny*cod.Nx + nx + 1] - p[ny*cod.Nx + nx - 1]);
				cod.v[ny*cod.Nx + nx] = dx_inv*(p[(ny + 1)*cod.Nx + nx] - p[(ny - 1)*cod.Nx + nx]);
			}
	}

	template<class T>
	void NS_Divergence(const NS_Vector &dom, ZNAC::LA::IVector<T> &cod, double dx)
	{
		const double dx_inv = 1/(2*dx);
		for(unsigned int ny = 1; ny < dom.Ny - 1; ++ny)
			for(unsigned int nx = 1; nx < dom.Nx - 1; ++nx)
				cod[ny*dom.Nx + nx] = dx_inv*((dom.u[ny*dom.Nx + nx + 1] - dom.u[ny*dom.Nx + nx - 1]) + (dom.v[(ny + 1)*dom.Nx + nx] - dom.v[(ny - 1)*dom.Nx + nx]));
	}
	template<class T>
	void NS_Laplace(const ZNAC::LA::IVector<T> &dom, ZNAC::LA::IVector<T> &cod, unsigned int Nx, unsigned int Ny, double dx)
	{
		const double dx2_inv = 1/(dx*dx);
		for(unsigned int ny = 1; ny < Ny - 1; ++ny)
			for(unsigned int nx = 1; nx < Nx - 1; ++nx)
				cod[ny*Nx + nx] = dx2_inv*(-4*dom[ny*Nx + nx] + ((dom[ny*Nx + nx + 1] + dom[ny*Nx + nx - 1]) + (dom[(ny + 1)*Nx + nx] + dom[(ny - 1)*Nx + nx])));
	}


	class NS_System
	{
	public:
		const unsigned int Nx, Ny;
		const double Re, a, min, max, dx;

		NS_Vector u;
		ZNAC::LA::Vector<double> p;

		NS_System(double Re, double a, unsigned int N, double width = 1);

		double x(int i);
		double y(int j);
	};

	//Staggered Cell
	//
	//    +------v_i,j+2-------+------v_i+1,j+2------+
	//    |                    |                     |
	//    |                    |                     |
	//  u_i,j+1  p_i,j+1  u_i+1,j+1   p_i+1,j+1 u\i+2,j+1
	//    |                    |                     |
	//    |                    |                     |
	//    +-------v_i,j+1------+------v_i+1,j+1------+
	//    |                    |                     |
	//    |                    |                     |
	//  u_i,j     p_i,j   u_i+1,j     p_i+1,j   u_i+2,j
	//    |                    |                     |
	//    |                    |                     |
	//    +-------v_i,j--------+------v_i+1,j--------+
	//

	struct NS_Index
	{
		const unsigned int Nx, Ny, margin;
		int &xi, &yi;

		NS_Index(int &xi, int &yi, unsigned int Nx, unsigned int Ny, unsigned int margin = 1):xi(xi), yi(yi), Nx(Nx), Ny(Ny), margin(margin){}
		unsigned int operator()(int dxi, int dyi) const {return (yi + dyi + margin)*(Nx + 2*margin) + xi + dxi + margin;}
		operator unsigned int()const{return (yi + margin)*(Nx + 2*margin) + xi + margin;}
	};

	class Diffusion
		:public ZNAC::LA::IMatrix<double>
	{
	public:
		Diffusion(double Dx, double Dy, unsigned int Nx, unsigned int Ny, unsigned int margin = 1);

		void operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const;
		double &operator()(unsigned int i, unsigned int j);
		const double &operator()(unsigned int i, unsigned int j)const;

	private:
		mutable int xi, yi;
		NS_Index n;
		const double Dx, Dy;
	};

	class NS_Grid
		:public ZNAC::LA::IVector<double>
	{
	public:
		const unsigned int Nx, Ny, margin;

		NS_Grid(unsigned int Nx, unsigned int Ny, unsigned int margin = 1);
		~NS_Grid();

		double &operator()(int xi, int yi);
		double &operator[](unsigned int i);
		const double &operator[](unsigned int i) const;
		unsigned int N()const;
		operator double*() const{return buf;}

	private:
		double *buf;
	};

	class StaggeredFunction
	{
	public:
		StaggeredFunction(double dx, double dy):dx_inv(1/dx), dy_inv(1/dy){}
	protected:
		const double dx_inv, dy_inv;
	};

	class Advection
		:public StaggeredFunction
	{
	public:
		Advection(double dx, double dy):StaggeredFunction(dx, dy){}
		void operator()(const NS_Grid &u, const NS_Grid &v, NS_Grid &du, NS_Grid &dv);
	};

	class Gradient
		:public StaggeredFunction
	{
	public:
		Gradient(double dx, double dy):StaggeredFunction(dx, dy){}
		void operator()(const NS_Grid &p, NS_Grid &dpu, NS_Grid &dpv);
	};

	class Divergence
		:public StaggeredFunction
	{
	public:
		Divergence(double dx, double dy):StaggeredFunction(dx, dy){}
		void operator()(const NS_Grid &u, const NS_Grid &v, NS_Grid &duv);
	};

	class NS_Coor
	{
	public:
		const double delta, origin;
		int i;

		NS_Coor(double delta, double origin);
		operator double();
		double operator()(int i);
	};

	class NS_Staggered
	{
	public:
		const double Re, dx, dy;
		NS_Grid u, v, p; 
		NS_Coor ux, uy, vx, vy, px, py;
		
		NS_Staggered(double Re, unsigned int Nx, unsigned int Ny, 
				double (*U)(double x, double y, int i, int j) = [](double x, double y, int i, int j)->double{return 0;}, 
				double (*V)(double x, double y, int i, int j) = [](double x, double y, int i, int j)->double{return 0;}, 
				double (*P)(double x, double y, int i, int j) = [](double x, double y, int i, int j)->double{return 0;}, 
				double left = -1, double right = 1, double bottom = -1, double top = 1, unsigned int margin = 1);

	};

	struct NS_Option
	{
		double Re;//Reynolds
		double a;//parabola parameter
		unsigned int N;//x-y divide num
		float width;//length from origin to bound
		double arrow_len;
		bool Advection;
		bool Diffusion;
		bool Force;
		bool Pressure;
		bool Auto;
		double dt;
		float wnd_width;
		NS_Option(int argc, char *argv[]);
		void print();
	};
}

#endif
