#ifndef __IPDE_1D__
#define __IPDE_1D__

#include <vector>
#include <string>
#include "CompressedM.hpp"

class PDE
{
public: // Instance vars

	enum Available_Solvers {
		Thomas,
		LU,
		LDU,
		Jacobi,
		Seidel,
		SOR,
		SSOR,
		GD,
		MR,
		CG
	};
 
public: //Functions

	virtual void Solve() = 0;
	virtual void Solve(const double omega) = 0;
	virtual void Solve_P() = 0;
	virtual void Solve_P(const double omega) = 0;
	virtual void ScaleSystem() = 0;
	virtual void ScaleSolution() = 0;
	virtual void Make_System() = 0;
	virtual void Save_Solution(const std::string &Filename) = 0;
	PDE(double a, double b, std::size_t N, const std::string &solve_method);
	PDE(double a, double b, std::size_t N, const std::string &solve_method, const std::string &Preconition_method);

protected: // Functions

	virtual double k(double x) = 0;
	virtual double f(double x) = 0;

protected: // Instance vars

	double a;
	double b;
	double h;
	std::size_t N;
	std::string Method;
	std::string Precond_Method;

	std::vector<double> Scaler;

	std::vector<double> u; // Solution
	std::vector<double> Rhs; // Right hand side of the system
	CMatrix Lhs; // Left hand side of the system
};

PDE::PDE(double a, double b, std::size_t N, const std::string &solve_method)
{
	this->a = a;
	this->b = b;
	this->N = N;
	this->Method = solve_method;
	h = (b - a) / (N - 1);
}

PDE::PDE(double a, double b, std::size_t N, const std::string &solve_method, const std::string &Precondition_method) : PDE(a,b,N,solve_method)
{
	this->Precond_Method = Precondition_method;
}

#endif