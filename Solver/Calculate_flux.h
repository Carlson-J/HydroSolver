//
// Created by carls502 on 1/20/2020.
//

#ifndef HYDROSOLVER_CALCULATE_FLUX_H
#define HYDROSOLVER_CALCULATE_FLUX_H

#include "System.h"
#include "Reconstruction.h"
#include "Riemann.h"
#include "Definitions.h"

// Tag for different ways to calculate flux
struct flux_linear_advection{};

// ------------------- Functor for calculating flux ------------- //
struct flux_functor{
    // Contains the grid
    const System system;

    explicit flux_functor(const System system) : system(system){};

    KOKKOS_INLINE_FUNCTION
    void operator() (const flux_linear_advection&, const long unsigned int& i) const{
        // We assume that this is a one dim grid with 1 guard cell.
        double left_state = system.grid(0,0, i-1, 0);
        double current_state = system.grid(0,0, i, 0);
        double dx = system.dx;
        double v = system.advection_velocity[DIRECTION_X];
        system.grid_flux(0, 0, 0, i, 0) = -v*(current_state - left_state)/dx;
    }

};

void Calculate_flux(const System& system){
    /*
     * To calculate fluxes we need to do a reconstruction and solve the Riemann problem
     */
    if (system.equation == "linear_advection" && system.time_stepping_algorithm == "forward_euler" && system.simulation_dim == 1){
        // For this we assume that there is 1 variable and one dim.
        Kokkos::RangePolicy<flux_linear_advection> policy(1,system.grid_flux.extent(X_DIM_FLUX)-1);
        Kokkos::parallel_for(policy, flux_functor(system));
    }
    else{
        Reconstruction(system);
        Riemann(system);
    }

}

#endif //HYDROSOLVER_CALCULATE_FLUX_H
