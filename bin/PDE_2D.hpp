#ifndef __IPDE_2D__
#define __IPDE_2D__

#include <vector>
#include <string>
#include "CompressedM.hpp"

class IPDE_2D
{
private:

protected: // Instance vars

	CMatrix Lhs;
	std::vector<double> Rhs;
	std::vector<double> u;

	std::string Solve_method;
	std::string Precondition_method;

	std::size_t N = (Nx) * (Ny);
	std::size_t Nx;
	std::size_t Ny;
	double hx;
	double hy;

public: // Functions

	virtual void MakeSystem() = 0;
	virtual void Solve(const double omega = 0) = 0;
	virtual void Solve_P(const double omega = 0) = 0;
	void SaveSolution(const std::string &);
	void SetSolver(const std::string &Solve_method);
	void SetPreconditioner(const std::string &Precondition_method);
	IPDE_2D(const std::size_t Nx, const std::size_t Ny, const std::string &Solve_Method);
	IPDE_2D(const std::size_t Nx, const std::size_t Ny, const std::string &Solve_method, 
			const std::string &Precondition_method);
	~IPDE_2D();

protected:

	virtual double f(const std::size_t k) const = 0;
	virtual double kxy(const std::size_t k) const = 0;
};

IPDE_2D::IPDE_2D(const std::size_t Nx, const std::size_t Ny, const std::string &Solver_Method)
{
	this->Nx = Nx;
	this->Ny = Ny;
	this->hx = (1 / (double)(Nx - 1));
	this->hy = (1 / (double)(Ny - 1));
	this->Solve_method = Solver_Method;
	this->N = Nx * Ny;
	Lhs.resize(N);
	Rhs.resize(N);
	u.resize(N);
}

IPDE_2D::IPDE_2D(const std::size_t Nx, const std::size_t Ny, const std::string &Solve_method, 
				 const std::string &Precondition_method) : IPDE_2D(Nx, Ny, Solve_method)
{
	this->Precondition_method = Precondition_method;
}

IPDE_2D::~IPDE_2D()
{
}

void IPDE_2D::SaveSolution(const std::string &Filename)
{
	std::ofstream out;
	out.open("gui/" + Filename + ".csv");
	if (out.is_open())
	{
		for (auto el : u)
		{
			out << el << std::endl;
		}
	}
}


#endif