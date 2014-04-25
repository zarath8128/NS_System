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
		void operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod);
	};

	class NS_Advection
		:public NS_Matrix
	{
	public:
		const double dx_inv;

		NS_Advection(unsigned int Nx, unsigned int Ny, double dx);
		double &operator()(unsigned int r, unsigned int c);
		const double &operator()(unsigned int r, unsigned int c)const;
		void operator()(const ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod);
	};

	void NS_Gradient(const ZNAC::LA::IVector<double> &p, NS_Vector &cod);
	void NS_Divergence(const NS_Vector &dom, ZNAC::LA::IVector<double> &cod);

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
}

#endif
