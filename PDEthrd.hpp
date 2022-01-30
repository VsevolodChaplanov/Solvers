#ifndef __PDEthrd_bc__
#define __PDEthrd_bc__

#include "PDE.hpp"
#include "CompressedM.hpp"
#include "Solvers.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

class PDEthrd : public PDE
{
public:

	PDEthrd(double a, double b, std::size_t N, const std::string &method);
	PDEthrd(double a, double b, std::size_t N, const std::string &method, double beta_left, double u_left, double beta_right, double u_right);
	void Solve() override;
	void Solve(const double omega) override;
	void Solve_P() override;
	void Solve_P(const double omega) override;
	void Save_Solution(const std::string &Filename) override;
	void SetSolver(const std::string &solve_method);
	void Make_System() override;
	void ScaleSystem() override;
	void ScaleSolution() override;
	void MakeExactSolution();
	void NormMaxRelExact();
	void CheckMatrixSymmetry();

protected:

	double f(double x) override;
	double k(double x) override;
	double u_ex(double x);

protected:

	double u_l = 0;
	double u_r = 1;
	double beta_l = 100000;
	double beta_r = 100000;
	// CMatrix Lhs;
	// std::vector<double> Rhs;
	// std::vector<double> u;
	std::vector<double> u_exact;
};


double PDEthrd::f(double x)
{
	return -16 * M_PI * M_PI * sin(16 * M_PI * x);
}

double PDEthrd::k(double x)
{
	return 1;
}

double PDEthrd::u_ex(double x)
{
	return x - sin(16 * M_PI * x) / 16 ;
}

PDEthrd::PDEthrd(double a, double b, std::size_t N, const std::string &method) : PDE(a,b,N, method) 
{ 
	Lhs.resize(N);
	u.resize(N);
	Rhs.resize(N);
}
PDEthrd::PDEthrd(double a, double b, std::size_t N, const std::string &method, double beta_left, double u_left, double beta_right, double u_right) : PDEthrd(a,b,N, method)
{
	this->u_l = u_left;
	this->u_r = u_right;
	this->beta_l = beta_left;
	this->beta_r = beta_right;
	Lhs.resize(N);
	u.resize(N);
	Rhs.resize(N);
}

void PDEthrd::Solve()
{
	IMatrixSolver* solver = IMatrixSolver::Fabric(Method);
	solver->solve(Lhs, Rhs, u);
	delete solver;
	solver = nullptr;
}

void PDEthrd::Make_System()
{
	Lhs.SetValue(0,0, k(h*1/2) / (h * h) + beta_l / h);
	Lhs.SetValue(0,1, - k(h*1/2) / (h * h));
	Lhs.SetValue(N-1,N-1, k(h*(N-1/2)) / (h * h) + beta_r / h);
	Lhs.SetValue(N-1,N-2, - k(h*(N-1/2)) / (h * h));

	Rhs[0] = f(0) / 2 + beta_l * u_l / h;
	Rhs[N-1] = f((N-1)*h) / 2 + beta_r * u_r / h;

	for (size_t i = 1; i < N-1; i++)
	{

		Lhs.SetValue(i,i, 2*k(i*h)/(h*h));
		Lhs.SetValue(i,i-1, -k(i*h)/(h*h));
		Lhs.SetValue(i,i+1, -k(i*h)/(h*h));

		Rhs[i] = f(i*h);	
	}
}

void PDEthrd::Save_Solution(const std::string &Filename)
{
	std::ofstream out;
	out.open(Filename + ".csv");
	if (out.is_open())
	{
		for (auto el : u)
		{
			out << el << std::endl;
		}
	}
}

void PDEthrd::Solve_P()
{
	if (Precond_Method.size() != 0)
	{
		IMatrixSolver* solver_p = IMatrixSolver::Fabric_P(Method, Precond_Method);
		solver_p->solve(Lhs, Rhs, u);
	}
	else
	{
		std::cout << "В конструкторе не был объявлен метод предобсулавливания, вызван не предобсуловленный метод решения" << std::endl;
		IMatrixSolver* solver = IMatrixSolver::Fabric(Method);
		solver->solve(Lhs, Rhs, u);
		//this->Solve();
	}
}

void PDEthrd::Solve(const double omega)
{
	IMatrixSolver* solver = IMatrixSolver::Fabric(Method, omega);
	solver->solve(Lhs, Rhs, u);
	delete solver;
	solver = nullptr;
}

void PDEthrd::Solve_P(const double omega)
{
	if (Precond_Method.size() != 0)
	{
		IMatrixSolver* solver_p = IMatrixSolver::Fabric_P(Method, Precond_Method);
		solver_p->solve(Lhs, Rhs, u);
	}
	else
	{
		std::cout << "В конструкторе не был объявлен метод предобсулавливания, вызван не предобсуловленный метод решения" << std::endl;
		IMatrixSolver* solver = IMatrixSolver::Fabric(Method, omega);
		solver->solve(Lhs, Rhs, u);
		//this->Solve();
	}
}

void PDEthrd::ScaleSystem()
{
	Scaler.resize(N);
	for (size_t i = 0; i < N; i++)
	{
		Scaler[i] = Lhs.GetValue(i,i);
	}
	
	for (size_t i = 0; i < N; i++)
	{
		Lhs.SetValue(i,i , 1);
		Rhs[i] = Rhs[i] / sqrt(Scaler[i]);
		for (auto elem : Lhs[i])
		{
			if (i != elem.first)
			{
				Lhs.SetValue(i, elem.first, elem.second / sqrt(Scaler[i]*Scaler[elem.first]));
			}
		}
	}
}

void PDEthrd::ScaleSolution()
{
	for (size_t i = 0; i < N; i++)
	{
		u[i] = u[i] / sqrt(Scaler[i]);
	}
}

void PDEthrd::MakeExactSolution()
{
	u_exact.resize(N);
	for (size_t i = 0; i < N; i++)
	{
		u_exact[i] = u_ex(h * i);
	}
}

void PDEthrd::NormMaxRelExact()
{
	double result = 0;
	for (size_t i = 0; i < N; i++)
	{
		double diff = fabs(u_exact[i] - u[i]);
		if (diff > result)
		{
			result = diff;
		}
	}
	std::cout << "Отклонение от точного решения: " << result << std::endl;
}

void PDEthrd::CheckMatrixSymmetry()
{
	if (CheckMatSym(Lhs) == false)
	{
		std::cout << "Матрица не симметричная" << std::endl;
	}
}

#endif