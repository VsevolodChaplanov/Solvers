#ifndef __IPDE_2D_3__
#define __IPDE_2D_3__

#include "PDE_2D.hpp"
#include "Solvers.hpp"
#include "Preconditioners.hpp"

class PDE_2D_third_bc : public IPDE_2D
{
protected:

	// CMatrix Lhs; // Матрица левой части
	// std::vector<double> Rhs; // Вектор правой части
	// std::vector<double> u; // Решение

	// Параметры метода решения
	// std::string Solve_method;
	// std::string Precondition_method;

	// Параметры сетки 
	// std::size_t Nx;
	// std::size_t Ny;
	// std::size_t N = Nx * Ny;
	// double hx = 1 / (Nx - 1);
	// double hy = 1 / (Ny - 1); 

	// Значения параметров в краевых условиях
	double u_low = 0;
	double beta_low = 100000;
	double u_up = 0;
	double beta_up = 100000;
	double u_right = 0;
	double beta_right = 100000;
	double u_left = 0;
	double beta_left = 100000;

	std::vector<double> Scaler;

public:

	void MakeSystem() override;
	void Solve(const double omega = 0) override;
	void Solve_P(const double omega = 0) override;
	void ScaleSystem();
	void ScaleSolution();
	PDE_2D_third_bc(const double beta_low, const double u_low, const double beta_left, const double u_left,
					const double beta_up, const double u_up, const double beta_right, const double u_right,
					const std::size_t Nx, const std::size_t Ny, const std::string &Solve_method, 
					const std::string &Precondition_method);
	PDE_2D_third_bc(const double beta_low, const double u_low, const double beta_left, const double u_left,
					const double beta_up, const double u_up, const double beta_right, const double u_right,
					const std::size_t Nx, const std::size_t Ny, const std::string &Solve_method);
	~PDE_2D_third_bc();

protected:

	double f(const std::size_t k) const override;
	double kxy(const std::size_t k) const override;
};

void PDE_2D_third_bc::Solve(const double omega)
{
	if (omega != 0)
	{
		IMatrixSolver* solver = IMatrixSolver::Fabric(Solve_method, omega);
		solver->solve(Lhs, Rhs, u);
		delete solver;
		solver = nullptr;
	}
	IMatrixSolver* solver = IMatrixSolver::Fabric(Solve_method);
	solver->solve(Lhs, Rhs, u);	
}

void PDE_2D_third_bc::Solve_P(const double omega)
{
	if (omega != 0 && Precondition_method.size() != 0)
	{
		IMatrixSolver* solver = IMatrixSolver::Fabric_P(Solve_method, Precondition_method, omega);
		solver->solve(Lhs, Rhs, u);
		delete solver;
		solver = nullptr;
	}	
	else if (Precondition_method.size() == 0)
	{
		std::cout << "Метод предобуславливания не был указан, используется не предобсуловленный метод решения" << std::endl;
		this->Solve(omega);
		return void();
	}
	IMatrixSolver* solver = IMatrixSolver::Fabric_P(Solve_method, Precondition_method);
	solver->solve(Lhs, Rhs, u);	
}

void PDE_2D_third_bc::MakeSystem()
{
	for (size_t k = 0; k < N; k++)
	{
		if (k == 0)
		{
			// Условия в точке (0,0)
			Lhs.SetValue(k, k, (2 / hx / hx * kxy(0) + 2 * beta_low / hx + 2 * beta_left / hy + 2 / hy / hy * kxy(0)) / 2);
			Lhs.SetValue(k, k+1, (- 2 / hx / hx * kxy(0))/2);
			Lhs.SetValue(k, k + (Nx), (- 2 / hy / hy * kxy(0))/2);
			Rhs[k] = (f(k) + 2 * beta_low / hx * u_low + 2 * beta_left / hy * u_left)/2;
		} else if (k == Nx - 1)
		{
			// Улсовия в точке (1, 0)
			Lhs.SetValue(k, k-1, (- 2 / hx / hx * kxy(0))/2);
			Lhs.SetValue(k, k, (2 / hx / hx * kxy(0) + 2 * beta_up / hx + 2 * beta_left / hy + 2 / hy / hy * kxy(0))/2);
			Lhs.SetValue(k, k + (Nx), (- 2 / hy / hy * kxy(0))/2);
			Rhs[k] = (f(k) + 2 * beta_low / hx * u_low + 2 * beta_right / hy * u_right)/2;
		} else if (k == N - (Nx))
		{
			// Условия в точке (0, 1)
			Lhs.SetValue(k, k-(Nx), (- 2 / hy / hy * kxy(0))/2);
			Lhs.SetValue(k, k, (2 / hx / hx * kxy(0) + 2 * beta_low / hx + 2 * beta_right / hy + 2 / hy / hy * kxy(0))/2);
			Lhs.SetValue(k, k+1, (- 2 / hx / hx * kxy(0))/2);
			Rhs[k] = (f(k) + 2 * beta_up / hx * u_up + 2 * beta_left / hy * u_left)/2;
		} else if (k == N - 1)
		{
			// Условия в точке (1, 1)
			Lhs.SetValue(k, k - (Nx), (- 2 / hy / hy * kxy(0)) / 2);
			Lhs.SetValue(k, k-1, (- 2 / hx / hx * kxy(0)) / 2);
			Lhs.SetValue(k, k, (2 / hx / hx * kxy(0) + 2 * beta_up / hx + 2 * beta_right / hy + 2 / hy / hy * kxy(0)) / 2);
			Rhs[k] = (f(k) + 2 * beta_up / hx * u_up + 2 * beta_right / hy * u_right) / 2;
		} else if (k < Nx && k != 0)
		{
			// Условия на нижней границе
			Lhs.SetValue(k, k-1, - 1 / hx / hx * kxy(0));
			Lhs.SetValue(k, k, 1 / hx / hx * kxy(0) + 1 / hx / hx * kxy(0) + 2 * beta_low / hy + 2 / hy / hy  * kxy(0));
			Lhs.SetValue(k, k+1, - 1 / hx / hx * kxy(0));
			Lhs.SetValue(k, k + (Nx), - 2 / hy /hy * kxy(0));
			Rhs[k] = f(k) + 2 * beta_low / hy * u_low;
		} else if (k > N - (Nx) && k != (N - 1))
		{
			// Верхняя граница
			Lhs.SetValue(k, k - (Nx), - 2 / hy /hy * kxy(0));
			Lhs.SetValue(k, k - 1, - 1 / hx / hx * kxy(0));
			Lhs.SetValue(k, k, 1 / hx / hx * kxy(0) + 1 / hx / hx * kxy(0) + 2 * beta_up / hy + 2 / hy / hy  * kxy(0));
			Lhs.SetValue(k, k + 1, - 1 / hx / hx * kxy(0));
			Rhs[k] = f(k) + 2 * beta_up / hy * u_up;
		} else if (k % (Nx) == 0 && k != 0 && k != N - (Nx))
		{
			// Левая граница
			Lhs.SetValue(k, k - (Nx), - 1 / hy /hy * kxy(0));
			Lhs.SetValue(k, k, 2 * beta_left / hx + 2 / hx / hx * kxy(0) + 1 / hy / hy * kxy(0) + 1 / hy / hy  * kxy(0));
			Lhs.SetValue(k, k + 1, - 2 / hx / hx * kxy(0));
			Lhs.SetValue(k, k + (Nx), - 1 / hy / hy * kxy(0));
			Rhs[k] = f(k) + 2 * beta_low / hx * u_left;
		} else if ((k + 1) % (Nx) == 0 && k != Nx && k != (N - 1))
		{
			// Правая граница
			Lhs.SetValue(k, k - (Nx), - 1 / hy /hy * kxy(0));
			Lhs.SetValue(k, k - 1, - 2 / hx / hx * kxy(0));
			Lhs.SetValue(k, k, 2 * beta_right / hx + 2 / hx / hx * kxy(0) + 1 / hy / hy * kxy(0) + 1 / hy / hy  * kxy(0));
			Lhs.SetValue(k, k + (Nx), - 1 / hy / hy * kxy(0));
			Rhs[k] = f(k) + 2 * beta_right / hx * u_right;
		} else
		{
			// Внутренняя область
			Lhs.SetValue(k, k - (Nx), 2 * (- 1 / hy / hy * kxy(0)));
			Lhs.SetValue(k, k - 1, 2 * (- 1 / hx / hx * kxy(0)));
			Lhs.SetValue(k, k, 2 * (1 / hx / hx * kxy(0) + 1 / hx / hx * kxy(0) + 1 / hy / hy * kxy(0) + 1 / hy / hy  * kxy(0)));
			Lhs.SetValue(k, k + 1, 2 * (-1 / hx / hx * kxy(0)));
			Lhs.SetValue(k, k + (Nx), 2 * (- 1 / hy / hy * kxy(0)));
			Rhs[k] = 2 * f(k);
		}
	}
}

double PDE_2D_third_bc::f(const std::size_t k) const
{
	
	int j = int(k / (Nx));
	int i = k - j * (Nx);
	double A = 5;
	double b = 1;
	double bx = 0.5, by = 0.5;
	return A * exp( - (( i*hx-bx ) * ( i*hx-bx ) + ( j*hy-by ) * ( j*hy-by )) * b );
}

double PDE_2D_third_bc::kxy(const std::size_t k) const
{
	return 1;
}

PDE_2D_third_bc::PDE_2D_third_bc(const double beta_low, const double u_low, const double beta_left, const double u_left,
								 const double beta_up, const double u_up, const double beta_right, const double u_right,
								 const std::size_t Nx, const std::size_t Ny, const std::string &Solve_method, 
								 const std::string &Precondition_method) : IPDE_2D(Nx, Ny, Solve_method, Precondition_method)
{ 
	this->beta_low = beta_low;
	this->beta_up = beta_up;
	this->beta_right = beta_right;
	this->beta_left = beta_left;
	this->u_low = u_low;
	this->u_left = u_left;
	this->u_up = u_up;
	this->u_right = u_right;
}

PDE_2D_third_bc::PDE_2D_third_bc(const double beta_low, const double u_low, const double beta_left, const double u_left,
								 const double beta_up, const double u_up, const double beta_right, const double u_right,
								 const std::size_t Nx, const std::size_t Ny, const std::string &Solve_method) : IPDE_2D(Nx, Ny, Solve_method)
{
	this->beta_low = beta_low;
	this->beta_up = beta_up;
	this->beta_right = beta_right;
	this->beta_left = beta_left;
	this->u_low = u_low;
	this->u_left = u_left;
	this->u_up = u_up;
	this->u_right = u_right;
}

void PDE_2D_third_bc::ScaleSystem()
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

void PDE_2D_third_bc::ScaleSolution()
{
	for (size_t i = 0; i < N; i++)
	{
		u[i] = u[i] / sqrt(Scaler[i]);
	}
}

PDE_2D_third_bc::~PDE_2D_third_bc()
{
}


#endif