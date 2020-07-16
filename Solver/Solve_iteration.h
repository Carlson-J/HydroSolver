//
// Created by carls502 on 1/20/2020.
//

#ifndef HYDROSOLVER_SOLVE_ITERATION_H
#define HYDROSOLVER_SOLVE_ITERATION_H

#include "System.h"
#include "Calculate_flux.h"
#include "Update_system.h"

void Solve_iteration(const System& system){
    /*
     * This will advance the system one time step.
     */
    // Choose which time stepping algorithm to use.
    if (system.time_stepping_algorithm == "forward_euler"){
        // This is a simple 1st order scheme, see IntroCompHydro eq. 4.2
        Calculate_flux(system);
        Update_system(system);
    }
    else if (system.time_stepping_algorithm == "direct"){
        /*
         * IntroCompHydro Eq 5.6-8 . We will reconstruct the once at the interfaces
         * half a time step forward in time.
         */
        Calculate_flux(system);
        Update_system(system);
    }
    else {
        throw(std::invalid_argument) {"Time stepping algorithm not found: " + system.time_stepping_algorithm};
    }
}

#endif //HYDROSOLVER_SOLVE_ITERATION_H
