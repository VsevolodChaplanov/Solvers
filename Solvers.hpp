#ifndef __MATRIX_SOLVERS__
#define __MATRIX_SOLVERS__

#include <vector>
#include <iostream>
#include "CompressedM.hpp"
#include "Preconditioners.hpp"
#include "VectorOperations.hpp"

/*-----------------------------Solvers Interface-----------------------------*/
class IMatrixSolver
{
protected:

	const std::size_t MAX_ITERATIONS = 10000000;
	const double eps = 0.0001;
	std::size_t Save_steps = 10;
	std::vector<double> R;

public:
	virtual void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) = 0;
	static IMatrixSolver* Fabric(const std::string&);
	static IMatrixSolver* Fabric(const std::string&, const double omega);
	static IMatrixSolver* Fabric_P(const std::string&, const std::string&);
	static IMatrixSolver* Fabric_P(const std::string&, const std::string&, const double omega);
};
/*-----------------------------Solvers Interface-----------------------------*/


/*-----------------------------Gradient decent solver-----------------------------*/
class GD_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Gradient decent solver-----------------------------*/


/*-----------------------------Minimal residuals solver-----------------------------*/
class MR_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Minimal residuals solver-----------------------------*/


/*-----------------------------Conjugate gradients solver-----------------------------*/
class CG_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Conjugate gradients solver-----------------------------*/


/*-----------------------------Symmetric Successive over relaxation solver-----------------------------*/
class SSOR_Solver : public IMatrixSolver
{
protected:

	double omega = 1.95;

public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
	void SetOmega(double omega);
};
/*-----------------------------Symmetric Successive over relaxation solver-----------------------------*/


/*-----------------------------Successive over relaxation solver-----------------------------*/
class SOR_Solver : public IMatrixSolver
{
protected:

	double omega = 1.95;

public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
	void SetOmega(double omega);
};
/*-----------------------------Successive over relaxation solver-----------------------------*/


/*-----------------------------Seidel-Gauss solver-----------------------------*/
class Seidel_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Seidel-Gauss solver-----------------------------*/


/*-----------------------------Jacobi solver-----------------------------*/
class Jacobi_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Jacobi solver-----------------------------*/


/*-----------------------------Thomas solver-----------------------------*/
class Thomas_Solver : public IMatrixSolver
{
protected:

	bool CheckTridiagonal(CMatrix& Lhs);

public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Thomas solver-----------------------------*/


/*-----------------------------Cholesky factorization solver-----------------------------*/
class LU_Solver : public IMatrixSolver
{
protected:

	void LU_decomposition(CMatrix &Lhs, CMatrix &L);

public:

	void solve(CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u) override;

};
/*-----------------------------Cholesky factorization solver-----------------------------*/


/*-----------------------------LDU factorization solver-----------------------------*/
class LDU_Solver : public IMatrixSolver
{
protected:

	void LDU_decomposition(CMatrix &Lhs, CMatrix &L, CMatrix &D);

public:

	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------LDU factorization solver-----------------------------*/



// Возможно есть смысл изначально унаследовать отдельный интерфейс для методов предобсуловленных

/*-----------------------------Gradient decent solver (P)-----------------------------*/
class GD_Solver_P : public IMatrixSolver
{
protected:
	IPreconditioner* Preconditioner;
public:

	GD_Solver_P(IPreconditioner* Preconditioner);
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Gradient decent solver (P)-----------------------------*/


/*-----------------------------Conjugate gradients solver (P)-----------------------------*/
class CG_Solver_P : public IMatrixSolver
{
protected:
	IPreconditioner* Preconditioner;
public:

	CG_Solver_P(IPreconditioner* Preconditioner);
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Conjugate gradients solver (P)-----------------------------*/


/*-----------------------------Realizations-----------------------------*/
IMatrixSolver* IMatrixSolver::Fabric(const std::string& method_name)
{
	if (method_name == "GD")
	{
		return new GD_Solver();
	} 
	else if (method_name == "MR")
	{
		return new MR_Solver();
	} 
	else if (method_name == "CG")
	{
		return new CG_Solver();
	} 
	else if (method_name == "SSOR")
	{
		return new SSOR_Solver();
	} 
	else if (method_name == "SOR")
	{
		return new SOR_Solver();
	} 
	else if (method_name == "Seidel")
	{
		return new Seidel_Solver();
	}
	else if (method_name == "Jacobi")
	{
		return new Jacobi_Solver();
	}
	else if (method_name == "LU")
	{
		return new LU_Solver();
	}
	else if (method_name == "LDU")
	{
		return new LDU_Solver();
	}
	return new Thomas_Solver();
}

IMatrixSolver* IMatrixSolver::Fabric(const std::string &method_name, const double omega)
{
	if (method_name == "SSOR")
	{
		SSOR_Solver* solver = new SSOR_Solver();
		solver->SetOmega(omega);
		return solver;
	} else if (method_name == "SOR")
	{
		SOR_Solver* solver = new SOR_Solver();
		solver->SetOmega(omega);
		return solver;
	}
	SOR_Solver* solver = new SOR_Solver();
	solver->SetOmega(omega);
	return solver;
}

IMatrixSolver* IMatrixSolver::Fabric_P(const std::string &method_name, const std::string &Precondition_method)
{
	IPreconditioner* Preconditioner = IPreconditioner::Fabric(Precondition_method);
	if (method_name == "CG")
	{
		return new CG_Solver_P(Preconditioner);
	}
	else if(method_name == "GD")
	{
		return new GD_Solver_P(Preconditioner);
	}
	return new GD_Solver_P(Preconditioner);
}

IMatrixSolver* IMatrixSolver::Fabric_P(const std::string &method_name, const std::string &Precondition_method, const double omega)
{
	IPreconditioner* Preconditioner = IPreconditioner::Fabric(Precondition_method, omega);
	if (method_name == "CG")
	{
		return new CG_Solver_P(Preconditioner);
	}
	else if(method_name == "GD")
	{
		return new GD_Solver_P(Preconditioner);
	}
	return new GD_Solver_P(Preconditioner);
}

void CG_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	double start_time = clock();
	//size_t N = Lhs.size();
	std::vector<double> r;
	std::vector<double> p;
	std::vector<double> R;

	r = VDiff((Lhs * u), Rhs);
	p = r;
	double Norm = DotProduct(r,p);
	std::vector<double> Ap = Lhs * p;
	double pAp = DotProduct(p, Ap);
	double lambda = Norm / pAp;
	u = VDiff(u, Mult_N(p, lambda));

	for (size_t k = 1; k < MAX_ITERATIONS; k++)
	{
		r = VDiff(r, Mult_N(Ap, lambda));
		double NewNorm = DotProduct(r,r);
		double alpha = NewNorm / Norm;
		p = VSum(r, Mult_N(p, alpha));
		Ap = Lhs * p;
		pAp = DotProduct(p, Ap);
		lambda = NewNorm / pAp;
		u = VDiff(u, Mult_N(p, lambda));
		Norm = NewNorm;
		if (k % Save_steps == 0)
		{
			R.push_back(sqrt(Norm));
		}
		if(sqrt(Norm) < eps * eps)
		{
			std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
			std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			std::cout << "Норма невязки: " << sqrt(Norm) << std::endl;
			WriteInFile(R, "Residuals");	
			break;
		}
	}
}

void GD_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	double start_time = clock(); 				// Замер времени
	std::vector<double> r;
	std::vector<double> R; 					// Для отслеживания нормы невязки с итерацией

	for (size_t k = 0; k < MAX_ITERATIONS; k++)
	{   
		r = VDiff(Lhs * u, Rhs);
		double Norm = DotProduct(r,r);
		double gamma = Norm / (DotProduct(Lhs * r, r));

		u = VDiff(u, Mult_N(r, gamma));

		if(k % Save_steps == 0) 
		{
			R.push_back(sqrt(Norm));
		} 
		
		// Условие выхода из цикла
		if(sqrt(Norm) < (eps))
		{
			std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
			std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << " секунд" << std::endl;
			WriteInFile(R, "Residuals");
			break;
		}
	} // k
}

void MR_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
    double start_time = clock();
    std::vector<double> r(Rhs.size());
    // Для отслеживания нормы невязки с итерацией
    std::vector <double> R;

    for (size_t k = 0; k < MAX_ITERATIONS; k++)
    {
        r = VDiff(Lhs * u, Rhs);

        std::vector<double> Ar = Lhs * r;
        double gamma = DotProduct(r , Ar) / DotProduct(Ar , Ar);
        u = VDiff(u, Mult_N(r, gamma));
        
        double Norm = sqrt(DotProduct(r,r));

        if(k % Save_steps == 0)
        {
            R.push_back(Norm);
        } 

        if(Norm < eps)
        {
            std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
            std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			WriteInFile(R, "Residuals");	
            break;
        }
    }
}

void SSOR_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	size_t N = Lhs.size();
	std::vector<double> Temp_U(N, 0.);
	std::vector<double> r(N, 0.);

	double start_time = clock();
	for (size_t k = 0; k < MAX_ITERATIONS; k++)
	{
		// Forward Propogation
		for (size_t i = 0; i < N; i++)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto el : Lhs[i])
			{
				if (el.first < i)
				{
					LowerSum += el.second * Temp_U[el.first];
				} // L * Temp_U
				if (el.first > i)
				{
					HigherSum += el.second * u[el.first];
				} // U * u
			} // End U * u & L * Temp_U // el
			
			Temp_U[i] = omega / Lhs.GetValue(i,i) * ( - LowerSum - HigherSum + Rhs[i]) + (1 - omega) * u[i];
		} // End i - filling solution vector // i

		// Back Propogation
		for (int i = N - 1; i >= 0; i--)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto el : Lhs[i])
			{
				if (el.first < i)
				{
					LowerSum += el.second * Temp_U[el.first];
				} // L * Temp_U
				if (el.first > i)
				{
					HigherSum += el.second * u[el.first];
				} // U * u
			} // End U * u & L * Temp_U // el
			
			u[i] = omega / Lhs.GetValue(i,i) * ( - LowerSum - HigherSum + Rhs[i]) + (1 - omega) * Temp_U[i];
		} // End i - filling solution vector // i

		double res = SecondNorm(VDiff(Lhs * u, Rhs));

		if (k % Save_steps == 0)
		{
			R.push_back(res);
		}

		if (res< eps)
		{
			std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
			std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << " секунд" << std::endl;
			WriteInFile(R, "Residuals");
			break;
		}		
	} // Iterations cycle // k
}

void SOR_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	size_t N = Lhs.size();
	std::vector<double> Temp_U(N, 0.);
	std::vector<double> r(N, 0.);

	double start_time = clock();
	for (size_t k = 0; k < MAX_ITERATIONS; k++)
	{
		for (size_t i = 0; i < N; i++)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto el : Lhs[i])
			{
				if (el.first < i)
				{
					LowerSum += el.second * Temp_U[el.first];
				} // L * Temp_U
				if (el.first > i)
				{
					HigherSum += el.second * u[el.first];
				} // U * u
			} // End U * u & L * Temp_U // el
			
			Temp_U[i] = omega / Lhs.GetValue(i,i) * ( - LowerSum - HigherSum + Rhs[i]) + (1 - omega) * u[i];
		} // End i - filling solution vector // i


		double res = SecondNorm(VDiff(Lhs * u, Rhs));

		if (k % Save_steps == 0)
		{
			R.push_back(res);
		}

		if (res < eps)
		{
			std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
			std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << " секунд" << std::endl;
			WriteInFile(R, "Residuals");
			break;
		}

		u = Temp_U;		
	} // Iterations cycle // k
}

void SSOR_Solver::SetOmega(double omega)
{
	this->omega = omega;
}

void SOR_Solver::SetOmega(double omega)
{
	this->omega = omega;
}

void Thomas_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	if (CheckTridiagonal(Lhs) == false)
	{
		std::cout << "Система не является трёхдиагональной" << std::endl;
	}
	double start_time = clock();
	size_t N = Lhs.size();
	// Массивы содержащие числа прямого хода метода прогонки
	std::vector<double> P(N, 0.);
	// Массивы содержащие числа прямого хода метода прогонки
	std::vector<double> Q(N, 0.);
	P[0] = -Lhs.GetValue(0,1)/(Lhs.GetValue(0,0));
	Q[0] = Rhs[0]/(Lhs.GetValue(0,0));
	for (size_t i = 1; i < N; i++)
	{
		P[i]=-Lhs.GetValue(i,i+1)/(Lhs.GetValue(i,i)+Lhs.GetValue(i,i-1)*P[i-1]);
		Q[i]=(-Lhs.GetValue(i,i-1)*Q[i-1]+Rhs[i])/(Lhs.GetValue(i,i)+Lhs.GetValue(i,i-1)*P[i-1]);
	}
	u[N-1] = (-Lhs.GetValue(N-1,N-2) * Q[N-2] + Rhs[N-1])/(Lhs.GetValue(N-1,N-1) + P[N-2] * Lhs.GetValue(N-1, N-2));
	for (size_t i = N-1; i > 0; i--)
	{
		u[i-1] = P[i-1] * (u[i]) + Q[i-1];
	}
	std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << " секунд" << std::endl;
}

bool Thomas_Solver::CheckTridiagonal(CMatrix& Lhs)
{
	bool check_result = true;
	for (size_t i = 0; i < Lhs.size(); i++)
	{
		Lhs[i].size() <= 3 ? check_result = true : check_result = false;
	}
	return check_result;
}

void Jacobi_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	double start_time = clock();
	std::size_t N = Lhs.size();
	std::vector<double> u_temp(N, 0.);

	for (size_t k = 0; k < MAX_ITERATIONS; k++)
	{
		for (size_t i = 0; i < N; i++)
		{
			double InnerSum = 0;
			for (auto elem : Lhs[i])
			{
				if (elem.first == i)
				{
					continue;
				}
				InnerSum += elem.second * u[elem.first];
			}

			u[i] = (-InnerSum + Rhs[i]) / Lhs.GetValue(i,i);
		}

		double res = SecondNorm(VDiff(Lhs * u, Rhs));

		if (k % Save_steps == 0)
		{
			R.push_back(res);
		}
		

		if (res < eps)
		{
			std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
			std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << " секунд" << std::endl;
			WriteInFile(R, "Residuals");
			break;
		}
	}
}

void Seidel_Solver::solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	size_t N = Lhs.size();
	std::vector <double> U(N, 0.);
	double start_time = clock();

	for (size_t k = 0; k < MAX_ITERATIONS; k++)
	{
		for (size_t i = 0; i < N; i++)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto elem : Lhs[i])
			{
				if (elem.first < i)
				{
				  	LowerSum += elem.second * U[elem.first];
				}
				if (elem.first > i)
				{
					HigherSum += elem.second * u[elem.first];
				}
			}
			
			U[i] = (- HigherSum - LowerSum + Rhs[i]) / Lhs.GetValue(i,i);
		}


		double res = SecondNorm(VDiff(Lhs * u, Rhs));

		if (k % Save_steps == 0)
		{
			R.push_back(res);
		}

		if (res < eps)
		{
			std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
			std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << " секунд" << std::endl;
			WriteInFile(R, "Residuals");
			break;
		}
		u = U;
	}
	
}

void LU_Solver::LU_decomposition(CMatrix &Lhs, CMatrix &L)
{
	std::size_t N = Lhs.size();
	double start_time = clock();

	// Ход алгоритма для нахождения L(ower)
	for (size_t i = 0; i < N; i++) 
	{
		for (size_t j = 0; j <= i; j++) 
		{
			double sum = 0;
			for (size_t k = 0; k < j; k++)
			{
				sum += L.GetValue(i,k) * L.GetValue(j,k);
			}
			if (i == j)
			{
				L.SetValue(i,i, sqrt(Lhs.GetValue(i,i) - sum));
			}
			else
			{
				double Diff = Lhs.GetValue(i,j) - sum;
				if (Diff != 0)
				{
					L.SetValue(i,j, (1.0 / L.GetValue(j,j) * (Diff)));
				}
			}
		}
	}
	std::cout << "Процедура LU разложения заняла: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}

void LU_Solver::solve(CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	double start_time = clock();
	size_t N = Lhs.size();
	CMatrix L(N);
	LU_decomposition(Lhs, L);
	double start_clock_Solution = clock();
	std::vector<double> y(N);

	for (size_t i = 1; i <= N; i++)
	{
		double Sum = 0;
		for (size_t j = 1; j <= i - 1; j++)
		{
			Sum += L.GetValue(i-1, j-1) * y[j-1];
			//Sum += L[i-1][j-1] * y[j-1];
		}
		y[i-1] = 1 / L.GetValue(i-1,i-1) * (Rhs[i-1] - Sum);
		//y[i-1] = 1 / L[i-1][i-1] * (RHS[i-1] - Sum);
	}

	for (size_t i = N; i > 0; i--)
	{
		double Sum = 0;
		for (size_t j = i+1; j < N + 1; j++)
		{
			Sum += L.GetValue(j-1,i-1) * u[j-1];
			//Sum += L[j-1][i-1] * x[j-1];
		}
		u[i-1] = 1 / L.GetValue(i-1,i-1) * (y[i-1] - Sum);
		//x[i-1] = 1 / L[i-1][i-1] * (y[i-1] - Sum);
	}

	std::cout << "Процедура решения СЛАУ с импользование LU разложения заняла: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
	std::cout << "Время решения без учета процедуры LU разложения: " << (clock() - start_clock_Solution) / CLOCKS_PER_SEC << std::endl;
}

void LDU_Solver::LDU_decomposition(CMatrix &Lhs, CMatrix &L, CMatrix &D)
{
	double start_time = clock();
	// Ход алгоритма для нахождения D(iagonal) и L(ower)
	for (size_t i = 0; i < Lhs.size(); i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			double sum = 0;
			for (size_t k = 0; k < j; k++)
			{
				sum += L.GetValue(i,k) * L.GetValue(j,k) * D.GetValue(k,k);
			}
			D.SetValue(i,i, Lhs.GetValue(i,i) - sum);
			L.SetValue(i,i,1);
			double diff = (Lhs.GetValue(i,j) - sum);
			if (diff != 0)
			{
				L.SetValue(i,j, 1 / D.GetValue(j,j) * diff);
			}
		}
	}
	std::cout << "Процедура LDU разложения заняла: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}

void LDU_Solver::solve(CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	double start_time = clock();
	size_t N = Lhs.size();
	CMatrix L(N);
	CMatrix D(N);

	LDU_decomposition(Lhs, L, D);

	std::vector<double> y(N);
	std::vector<double> x(N);

	for (size_t i = 1; i <= N; i++)
	{
		double Sum = 0;
		for (size_t j = 1; j <= i - 1; j++)
		{
			Sum += D.GetValue(j-1,j-1) * L.GetValue(i-1,j-1) * y[j-1];
		}
		y[i-1] = 1 / D.GetValue(i-1,i-1) * (Rhs[i-1] - Sum);
	}

	for (size_t i = N; i > 0; i--)
	{
		double Sum = 0;
		for (size_t j = i+1; j < N + 1; j++)
		{
			Sum += L.GetValue(j-1,i-1) * u[j-1];
		}
		u[i-1] = 1 / L.GetValue(i-1, i-1) * (y[i-1] - Sum);
	}

	std::cout << "Процедура решения СЛАУ с импользование LDU разложения заняла: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}

void CG_Solver_P::solve(CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	double start_time = clock();
	size_t N = Lhs.size();
	std::vector<double> w;
	std::vector<double> r = VDiff(Lhs * u, Rhs);
	std::vector<double> p = r;
	std::vector<double> Ap = Lhs * p;
	double rw = DotProduct(r,r);
	double lambda = (DotProduct(p,r)) / (DotProduct(p, Ap));
	u = VDiff(u, Mult_N(p, lambda));

	for (size_t k = 0; k < MAX_ITERATIONS; k++)
	{	
		r = VDiff(r, Mult_N(Ap, lambda));
		w = Preconditioner->Precondition(Lhs, r);
		double RW = DotProduct(r, w);
		double Norm2 = DotProduct(r, r);
		double alpha = RW * (1 / rw);
		rw = RW;
		p = VSum(w, Mult_N(p, alpha));
		Ap = Lhs * p;
		double pAp = DotProduct(p, Ap);
		lambda = rw / pAp;
		u = VDiff(u, Mult_N(p, lambda));

        if(k % Save_steps == 0)
        {
            //Y_iter.push_back(Y);
            R.push_back(sqrt(Norm2));
        } 

		if(sqrt(Norm2) < eps * eps)
        {
            std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
            std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			WriteInFile(R, "Residuals_P");
            break;
        }
	}
}

void GD_Solver_P::solve(CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	//! Видимо в методичке ошибка
	double start_time = clock();
    std::vector<double> r;
	std::vector<double> w;
    // Для отслеживания нормы невязки с итерацией
    std::vector <double> R;

    for (size_t k = 0; k < MAX_ITERATIONS; k++)
    {   
        // r = Ay-B
        r = VDiff(Lhs * u, Rhs);
		w = Preconditioner->Precondition(Lhs, r);
        // gamma = (r,r) / (Ar,r)
		double Norm2 = DotProduct(r,r);
        double gamma = DotProduct(r, w) / DotProduct(r, Lhs * w);

        u = VDiff(u, Mult_N(r, gamma));

        if(k % Save_steps == 0)
        {
            //Y_iter.push_back(Y);
            R.push_back(sqrt(Norm2));
        } 

        if(sqrt(Norm2) < eps)
        {
            std::cout << "Процесс сошелся на итерации: " << k <<std::endl;
            std::cout << "Время затраченное на решение: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			WriteInFile(R, "Residuals_P");
            break;
        }
    }
}

CG_Solver_P::CG_Solver_P(IPreconditioner* Preconditioner)
{
	this->Preconditioner = Preconditioner;
}

GD_Solver_P::GD_Solver_P(IPreconditioner* Preconditioner)
{
	this->Preconditioner = Preconditioner;
}

#endif
