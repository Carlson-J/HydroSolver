////
//// Created by carls502 on 1/23/2020.
////
///*
// * This file contains the functions that can be used to benchmark the HydroSolver
// */
//
//#ifndef HYDROSOLVER_BENCHMARK_SIMULATIONS_H
//#define HYDROSOLVER_BENCHMARK_SIMULATIONS_H
//
//#include "string"
//#include "Solver/Timer.h"
//#include "Solver/Solver.h"
//#include <iomanip>
//
//std::string Linear_advection_benchmark(){
//    std::ostringstream output;
//    output << "Linear Advection Benchmark" << endl;
//    // Benchmarking parameters
//    int numRuns = 1;
//    double runTime = 1;
//    // System parameters
//    double dx = 0.001;
//    double cfl = 0.8;
//    double dt = 0.1;
//    int size_x = 10000;
//    int size_y = 1;
//    int size_z = 1;
//    int num_vars = 1;
//    string equation = "linear_advection";
//    string dt_type = "cfl_linear_advection";
//    string grid_initial = "advect_1d_hat";
//    string boundary_condition_x = "periodic";
//    string boundary_condition_y = "periodic";
//    string boundary_condition_z = "periodic";
//    string time_stepping_algorithm = "forward_euler";
//    double linear_advection_velocity = 1;
//    int save_frequency = 1;
//    bool intermediate_saving = true;
//    bool print_output = false;
//    for (int j = 10; j <= 1000; j = j*10){
//        output << "Number of cells: " << j << endl;
//        size_x = j;
//        // Create System
//        const System system(dx, cfl,dt, size_x , size_y, size_z, num_vars, dt_type, grid_initial,
//                      boundary_condition_x, boundary_condition_y, boundary_condition_z,
//                      time_stepping_algorithm, equation, linear_advection_velocity, save_frequency,
//                      intermediate_saving, print_output);
//        // Create Simulation
//        Solver solver(system);
//        // Run many times
//        for (int i = 0; i < numRuns; ++i){
//        // Time how long it takes to run
//            solver.timer.Start_timer("Linear Advection");
//            solver.Run_simulation(runTime);
//            solver.timer.Stop_time("Linear Advection");
//            output << solver.timer.Get_time("Linear Advection") << endl;
//        }
//
//    }
//    return output.str();
//}
//
//
//#endif //HYDROSOLVER_BENCHMARK_SIMULATIONS_H
