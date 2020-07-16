//
// Created by carls502 on 1/17/2020.
//
// The system structure holds all the information needed to run a simulation, including the grid.

#ifndef HYDROSOLVER_SYSTEM_H
#define HYDROSOLVER_SYSTEM_H

#include "../../kokkos/core/src/Kokkos_Core.hpp"
#include "Definitions.h"
using std::string;

struct System{
    // Debug settings
    bool print_output;

    // Equation to solve
    const string equation;

    // Grid
    int guardCells;
    Kokkos::View<double ****> grid;
    Kokkos::View<double *****> grid_flux_L; // The left side of the interface at i+1/2 on the grid will be the i value of this grid
    Kokkos::View<double *****> grid_flux_R; // The right side of the interface at i+1/2 on the grid will be the i+1 value of this grid
    Kokkos::View<double *****> grid_flux;
    string grid_initial_state;
    double dx; // Cell size in x direction
    int simulation_dim;
    // Velocities
    double advection_velocity[3];
    const string boundary_condition_x;
    const string boundary_condition_y;
    const string boundary_condition_z;

    // Time parameters and state variables
    Kokkos::View<double*> dt;
    double t;
    const double cfl;
    const string calc_dt_type;

    // Time stepping algorithm
    const string time_stepping_algorithm;

    // Reconstruction types
    const string reoncstruct_type = "direct_1D_taylor";

    // Riemann problem
    const string riemann_type = "simple_upwinding";

    // Limiters
    const string limiter;

    // Saving parameters
    int steps;
    int last_save_step;
    const int save_frequency;
    const bool intermediate_saving;

    System(double dx, double cfl, double dt_new, int sizeX, int sizeY, int sizeZ, int numVars,
            string calc_dt_type, string grid_initial_state, string boundary_condition_x,
            string boundary_condition_y, string boundary_condition_z, string time_stepping_algorithm,
            string equation, double linear_advection_velocity, int save_frequency,
            bool intermediate_saving, bool print_output, string limiter) :
        steps(0), last_save_step(0), t(0.0),
        cfl(cfl), dx(dx), dt("dt", 1),
        calc_dt_type(std::move(calc_dt_type)), grid_initial_state(std::move(grid_initial_state)),
        boundary_condition_x(std::move(boundary_condition_x)),
        boundary_condition_y(std::move(boundary_condition_y)),
        boundary_condition_z(std::move(boundary_condition_z)),
        time_stepping_algorithm(std::move(time_stepping_algorithm)), equation(equation),
        save_frequency(save_frequency),
        intermediate_saving(intermediate_saving), print_output(print_output),
        limiter(limiter){
        // Note that the grid is initialized here
        if (print_output) {std::cout << "System Initialized" << std::endl;}

        // Check to see if input is reasonable
        assert(cfl <= 1);
        assert(sizeX > 0 && sizeY > 0 && sizeZ > 0 && numVars > 0);
        assert(dx > 0);

        // Determine simulation dimension
        simulation_dim = 1;

        // Set number of guard cells and othre parametersbased on parameters
        if (equation == "linear_advection") {
            if (!(sizeX > 1 && sizeY == 1 & sizeZ == 1 && numVars == 1)) {
                throw (std::invalid_argument) {"For linear advection you must have a 1D grid with 1 variable"};
            }
            advection_velocity[DIRECTION_X] = linear_advection_velocity;
            advection_velocity[DIRECTION_Y] = 0;
            advection_velocity[DIRECTION_Z] = 0;
            if (this->time_stepping_algorithm == "direct") {
                guardCells = 2;
            } else if (this->time_stepping_algorithm == "forward_euler"){
                guardCells = 1;
            } else {
                throw (std::invalid_argument) {"Invalid time stepping algorithm used for linear advection problem."};
        }

        }
        else {
            guardCells = 0;
            advection_velocity[DIRECTION_X] = 0;
            advection_velocity[DIRECTION_Y] = 0;
            advection_velocity[DIRECTION_Z] = 0;
        }
        // Initialize grid
        // --- deterime dimensions based on guard cells needed and simulation dimension
        int Nx = sizeX + guardCells*2;
        int Ny = 1;
        int Nz = 1;
        if (simulation_dim > 1){
            Ny = sizeY + guardCells*2;
        }
        if (simulation_dim > 2){
            Nz = sizeZ + guardCells*2;
        }


        grid = Kokkos::View<double ****>("Computational Grid", Nz, Ny, Nx, numVars);
        grid_flux_L = Kokkos::View<double *****>("Flux_L Grid", simulation_dim, Nz, Ny, Nx, numVars);
        grid_flux_R = Kokkos::View<double *****>("Flux_R Grid", simulation_dim, Nz, Ny, Nx, numVars);
        grid_flux = Kokkos::View<double *****>("Flux Grid", simulation_dim, Nz, Ny, Nx, numVars);
    }

//    int Get_range
};


#endif //HYDROSOLVER_SYSTEM_H
