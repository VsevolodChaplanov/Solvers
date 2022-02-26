cmake_minimum_required(VERSION 3.8)

project(Solvers VERSION 1.0)

set(SOURCE_EXE bin/main.cpp)
set(HEADER bin/CompressedM.hpp
			bin/PDE_2D_third_bc.hpp
			bin/PDE_2D.hpp
			bin/PDE.hpp
			bin/PDEthrd.hpp
			bin/Preconditioners.hpp
			bin/Solvers.hpp
			bin/VectorOperations.hpp)

add_executable(main ${SOURCE_EXE})