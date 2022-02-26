#include "PDEthrd.hpp"
#include "PDE_2D_third_bc.hpp"
#include <string>
#include <memory>
#include <random>

// testline
int main(int argc, char const *argv[])
{
	// Одномерная задача
	// std::unique_ptr<PDEthrd> My_PDE(new PDEthrd(0,1,51, "SSOR"));
	// My_PDE->Make_System();
	// My_PDE->Solve();
	// My_PDE->Save_Solution("Solution");

	// Двумерная задача
	std::unique_ptr<PDE_2D_third_bc> My_PDE(new PDE_2D_third_bc(10000, 0, 10000, 0, 10000, 0, 10000, 0, 11, 11, "GD"));
	My_PDE->MakeSystem();
	//My_PDE->ScaleSystem();
	My_PDE->Solve();
	//My_PDE->ScaleSolution();
	My_PDE->SaveSolution("2D");

	// Случайная матрица
	// std::size_t N = 10000;
	// CMatrix RandomLhs(N, "Random");
	// std::vector<double> Solution(N);
	// for (size_t i = 0; i < N; i++)
	// {
	// 	Solution[i] = rand() % 5;
	// }

	// std::vector<double> Rhs = RandomLhs * Solution;
	// std::vector<double> N_Solution(N);

	// IMatrixSolver* solver = IMatrixSolver::Fabric_P("CG", "ISO0_2");
	// solver->solve(RandomLhs, Rhs, N_Solution);
	// std::cout << "Отклонение от модельного решения: " << NormMax(Solution, N_Solution) << std::endl;


	// // std::vector<double> N_Solution_cg(N);
	// // IMatrixSolver* solver_cg = IMatrixSolver::Fabric_P("CG", "ISO0_2");
	// // solver_cg->solve(RandomLhs, Rhs, N_Solution_cg);

	// // std::cout << "Отклонение от модельного решения: " << NormMax(Solution, N_Solution_cg) << std::endl;
	
	return 0;
}
