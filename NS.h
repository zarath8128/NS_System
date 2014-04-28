#ifndef ZARATH_NS_NS_H
#define ZARATH_NS_NS_H

#include <ZNAC/LA/Vector.h>
#include <ZNAC/LA/Matrix.h>

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
