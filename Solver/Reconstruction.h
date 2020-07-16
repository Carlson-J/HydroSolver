//
// Created by carls502 on 1/20/2020.
//

#ifndef HYDROSOLVER_RECONSTRUCTION_H
#define HYDROSOLVER_RECONSTRUCTION_H

#include "System.h"
#include "Iteration_space.h"
#include <cmath>

// Tags for different ways of calculating reconstruction
struct recon_1D_taylor_tag {};
struct recon_1D_taylor_minmod_tag{};


// ------------ Functor for calculating reconstruction given the system's state -------- //
struct Reconstruction_functor{

    const System system;

    explicit Reconstruction_functor(const System& system): system(system) {
    };

    KOKKOS_INLINE_FUNCTION
    void operator() (const recon_1D_taylor_tag&, const long unsigned int& i,const long unsigned int& v) const{
        // Reconstruct the i+1/2 cell interface
        // It is assumed that this is a 1D problem.
        // thus we set all other dim coordinates to the first entry.
        //IntoCompHydro eq. 5.6-5.8
        int g = system.guardCells;
        // Eq 5.6
        system.grid_flux_L(DIRECTION_X,0,0,i,v) = system.grid(0,0, i, v) +
                0.25 * (system.grid(0,0,i+1,v) - system.grid(0,0,i-1,v))
                * (1 - system.dt(0)/system.dx * system.advection_velocity[DIRECTION_X]);
        // Eq 5.7
        system.grid_flux_R(DIRECTION_X,0,0,i,v) = system.grid(0,0, i+1, v) -
                0.25 * (system.grid(0,0,i+2,v) - system.grid(0,0,i,v))
                * (1 + system.dt(0)/system.dx * system.advection_velocity[DIRECTION_X]);
    }
    KOKKOS_INLINE_FUNCTION
    void operator() (const recon_1D_taylor_minmod_tag&, const long unsigned int& i,const long unsigned int& v) const{
        // Reconstruct the i+1/2 cell interface
        // It is assumed that this is a 1D problem.
        // thus we set all other dim coordinates to the first entry.
        //IntoCompHydro eq. 5.6-5.8 & eq.5.11-5.12
        int g = system.guardCells;
        // Calculate temp dadx's (eq 5.11)
        double dadx_i = system.grid(0,0,i,v) - system.grid(0,0,i-1,v);
        double dadx_ip1 = system.grid(0,0,i+1,v) - system.grid(0,0,i,v);
        double dadx_ip2 = system.grid(0,0,i+2,v) - system.grid(0,0,i+1,v);
        // Determine dadx_i to use for i and i+1
        // If a x b < 0 then make the slope zero
        if(dadx_i*dadx_ip1 <= 0){
            dadx_i = 0;
        }
        else {
            // if |b| < |a| we use b (dadx_ip1) otherwise we use a (dadx_i so we keep it the same)
            if(std::fabs(dadx_ip1) < std::fabs(dadx_i)){
                dadx_i = dadx_ip1;
            }
            // if |b| < |a| we use b (dadx_ip2) otherwise we use a (dadx_ip1 so we keep it the same)
            if (std::fabs(dadx_ip2) < std::fabs(dadx_ip1)){
                dadx_ip1 = dadx_ip2;
            }
        }

        // Eq 5.6
        system.grid_flux_L(DIRECTION_X,0,0,i,v) = system.grid(0,0, i, v) +
                0.5 * dadx_i
                * (1 - system.dt(0)/system.dx * system.advection_velocity[DIRECTION_X]);
        // Eq 5.7
        system.grid_flux_R(DIRECTION_X,0,0,i,v) = system.grid(0,0, i+1, v) -
                0.5 * dadx_ip1
                * (1 + system.dt(0)/system.dx * system.advection_velocity[DIRECTION_X]);
    }
};


// ------------------ Callable Function for calculating time step ---------------------- //

void Reconstruction(const System& system){
    /*
     * Reconstruct system based on type of reconstruction specified
     * TODO: Pass a stencil to be reconstructed over
     */
    // Find starting and ending index based on simulation dimension
    Iteration_space itrSpc(system);
    // Create Functor
    const Reconstruction_functor recon_functor(system);
    // Choose reconstruction type
    if (system.reoncstruct_type == "direct_1D_taylor"){
        if (system.limiter == "minmod"){
            Kokkos::MDRangePolicy<recon_1D_taylor_minmod_tag, Kokkos::Rank<2>> policy(
                    {itrSpc.N_flux_x_0 , itrSpc.N_v_0},{itrSpc.N_flux_x_f, itrSpc.N_v_f}); // The +1 and -1 come from needing interfaces out outer cells
            // Use Kokkos and functor to initialize grid
            Kokkos::parallel_for(policy, recon_functor);
        }
        else{
            Kokkos::MDRangePolicy<recon_1D_taylor_tag, Kokkos::Rank<2>> policy(
                    {itrSpc.N_flux_x_0 , itrSpc.N_v_0},{itrSpc.N_flux_x_f, itrSpc.N_v_f}); // The +1 and -1 come from needing interfaces out outer cells
            // Use Kokkos and functor to initialize grid
            Kokkos::parallel_for(policy, recon_functor);
        }
    }

}

#endif //HYDROSOLVER_RECONSTRUCTION_H
