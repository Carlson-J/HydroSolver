//
// Created by carls502 on 1/17/2020.
//

#ifndef HYDROSOLVER_CALCULATE_DT_H
#define HYDROSOLVER_CALCULATE_DT_H

#include "System.h"


// Tag for different ways to calculate time step
struct calc_dt_cfl_tag{};


// ------------------- Functor for calculating dt given the system's state ------------- //
struct Calc_dt_functor{

    KOKKOS_INLINE_FUNCTION
    void operator() (const calc_dt_cfl_tag&, const int& i, double& dt) const{
        // Sets the time step based on convergence condition by cfl (Courant–Friedrichs–Lewy)
        //TODO: implement and reference clf
        if (i == 0 && dt < 1.0){
            dt += 1.0;
        }
    }

};

// ------------------ Callable Function for calculating time step ---------------------- //
double Calc_dt(double runTime, double time, const System& system){
    /*
     * This method is called each time step to calculate the size of the next time step.
     * The time stepping function used is set by Set_dt
     */
    // Choose which type of dt calculator to use
    double new_dt = 0;
    if (system.calc_dt_type == "cfl_linear_advection"){
        /* Linear advection has a constant time step based off of the velocity, which is constant for all cells, and
         * the cell size. See IntroCompHydro 4.1
         */
        new_dt = system.cfl * system.dx / fabs(system.advection_velocity[DIRECTION_X]);
    }
    else if (system.calc_dt_type == "fixed"){
        // Just keep the time step constant
        Kokkos::parallel_reduce(1, KOKKOS_LAMBDA (const int i, double& dt_tmp ){
            dt_tmp = system.dt(0);
        }, new_dt);
    }
    else if ( system.calc_dt_type == "cfl"){
        std::cerr << "clf not implemented. Setting dt to 0.5 of current dt." << std::endl;
        Kokkos::parallel_reduce(Kokkos::RangePolicy<calc_dt_cfl_tag>(0,1), Calc_dt_functor(), new_dt);
    }
    else {
        throw(std::invalid_argument){"Invalid dt_type used."};
    }

    // Check if the timestep goes past runtime. If so, fix it so it finished on time
    if (runTime < time + new_dt){
        new_dt = runTime - time;
    }
    // Save new dt in system
    Kokkos::parallel_for(1, KOKKOS_LAMBDA (const int i){
         system.dt(0) = new_dt;
    });
    return new_dt;
}






#endif //HYDROSOLVER_CALCULATE_DT_H
