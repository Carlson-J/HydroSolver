cmake_minimum_required(VERSION 3.10)
project(HydroSolver)

add_subdirectory(../kokkos kokkos_build)
include_directories(${Kokkos_INCLUDE_DIRS_RET})
include_directories(../include/json)

set(CMAKE_CXX_STANDARD 14)

add_executable(HydroSolver main.cpp Solver/Solver.h Solver/System.h Solver/Calculate_dt.h Solver/Initialize_grid.h Solver/boundary_condtions.h Solver/Definitions.h Solver/Solve_iteration.h Solver/Calculate_flux.h Solver/Reconstruction.h Solver/Riemann.h Solver/Update_system.h Solver/Timer.h Benchmark_simulations.h Test_simulations.h Solver/Iteration_space.h)
target_link_libraries(HydroSolver kokkos)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

# make compiled code possition independent. This is needed to wrap it in python (https://stackoverflow.com/questions/38296756/what-is-the-idiomatic-way-in-cmake-to-add-the-fpic-compiler-option)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
