//
// Created by carls502 on 1/20/2020.
//

#ifndef HYDROSOLVER_UPDATE_SYSTEM_H
#define HYDROSOLVER_UPDATE_SYSTEM_H

#include "System.h"
#include "Definitions.h"


void Update_system(const System& system){
    /*
     * Update the system based on fluxes
     *
     */
    // Find starting and ending index based on simulation dimension
    Iteration_space itrSpc(system);
    // Update each directions
    Kokkos::MDRangePolicy<Kokkos::Rank<4>> policy({itrSpc.N_z_0, itrSpc.N_y_0, itrSpc.N_x_0, itrSpc.N_v_0},
            {itrSpc.N_z_f, itrSpc.N_y_f, itrSpc.N_x_f, itrSpc.N_v_f});

    int offsets[3] = {0,0,0};
    for (int direction = 0; direction < system.simulation_dim; direction++){
        offsets[direction] = 1;
        if (system.time_stepping_algorithm == "forward_euler"){
            Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const int& k, const int& j, const int& i, const int& v){
                system.grid(k,j,i,v) += system.dt(0) * system.grid_flux(direction,k,j,i,v);
            });
        }
        else if (system.time_stepping_algorithm == "direct"){
            // IntroCompHydro Eq 5.4
            // Note that the ith flux corresponds to the i-1/2 interface.
            Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const int& k, const int& j, const int& i, const int& v){
                system.grid(k,j,i,v) -= system.dt(0)/system.dx *
                        (system.grid_flux(direction,k+offsets[2],j+offsets[1],i+offsets[0],v) - system.grid_flux(direction,k,j,i,v));
            });
        }
        else{
            throw(std::invalid_argument){"Invalid time_stepping_algorithm"};
        }
        offsets[direction] = 0;
        Kokkos::fence();
    }



}
#endif //HYDROSOLVER_UPDATE_SYSTEM_H
