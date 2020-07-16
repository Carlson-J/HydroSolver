//
// Created by carls502 on 1/20/2020.
//

#ifndef HYDROSOLVER_RIEMANN_H
#define HYDROSOLVER_RIEMANN_H

#include "System.h"
#include "Iteration_space.h"

// Tags for different ways of calculating reconstruction
struct riemann_upwinding_1D_tag {};


// ------------ Functor for calculating reconstruction given the system's state -------- //
struct Riemann_functor{

    const System system;

    explicit Riemann_functor(const System& system): system(system) {
    };

    KOKKOS_INLINE_FUNCTION
    void operator() (const riemann_upwinding_1D_tag&, const long unsigned int& i,const long unsigned int& v) const{
        // Use simple upwinding to update the flux grid using the left and right flux values. 1D only!
        // Note that the ith flux corresponds to the i-1/2 interface.
        //IntoCompHydro eq. 5.10
        int g = system.guardCells;
        if (system.advection_velocity[DIRECTION_X] > 0){
            // Right side of interface
            system.grid_flux(DIRECTION_X,0,0,i,v) = system.grid_flux_L(DIRECTION_X,0,0,i-1,v);
        } else{
            // Right side of interface
            system.grid_flux(DIRECTION_X,0,0,i,v) = system.grid_flux_R(DIRECTION_X,0,0,i,v);
        }
    }
};


void Riemann(const System& system){
    /*
     * Solves the Riemann problem and stores the results in the grid_flux view.
     */
    // Find starting and ending index based on simulation dimension
    Iteration_space itrSpc(system);
    // Create Functor
    const Riemann_functor riemann_functor(system);
    // Choose type of riemann solver to use
    if (system.riemann_type == "simple_upwinding" && system.simulation_dim == 1){
        Kokkos::MDRangePolicy<riemann_upwinding_1D_tag, Kokkos::Rank<2>> policy(
                {itrSpc.N_x_0, itrSpc.N_v_0},
                {itrSpc.N_x_f + 1, itrSpc.N_v_f}); // Note that there are +1 interfaces compared to cells (not including host cells)
        // Use Kokkos and functor to initialize grid
        Kokkos::parallel_for(policy, riemann_functor);
    } else{
        throw(std::invalid_argument){"Unable to find Riemann solver for given configuration."};
    }
}

#endif //HYDROSOLVER_RIEMANN_H
