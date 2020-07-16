//
// Created by carls502 on 1/17/2020.
//
#include "System.h"
#ifndef HYDROSOLVER_INITIALIZE_GRID_H
#define HYDROSOLVER_INITIALIZE_GRID_H

#include "Definitions.h"


// Tag for different ways to initialize grid
struct grid_initialize_testing_tag{};
struct grid_initialize_1d_advect_tag{};


// ------------------- Functor for calculating dt given the system's state ------------- //
struct Initialize_grid_functor{
    // Contains the grid
    const System system;

    explicit Initialize_grid_functor(const System& system) : system(system){};

    KOKKOS_INLINE_FUNCTION
    void operator() (const grid_initialize_testing_tag&, const long unsigned int& k, const long unsigned int& j, const long unsigned int& i,const long unsigned int& v) const{
        // Initialize grid to a state that is easy to examine.
        system.grid(k,j,i,v) = i + j;
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (const grid_initialize_1d_advect_tag&, const long unsigned int& k, const long unsigned int& j, const long unsigned int& i,const long unsigned int& v) const{
        // Initialize to a hat wave. See IntroCompHydro eq. 4.6
        double value = 0;
        if(i >= (1.0/3.0) * (system.grid.extent(X_DIM_GRID) - system.guardCells) && i <= (2.0/3.0) * system.grid.extent(X_DIM_GRID) ){
            value = 1;
        }
        system.grid(k, j, i, v) = value;
    }

};

// ------------------ Callable Function for initializing grid ---------------------- //
void Initialize_grid(const System& system){
    /*
     * Initialize grid to starting values
     */
    // Get grid dimensions
    int N_v = system.grid.extent(VAR_DIM_GRID);
    int N_x = system.grid.extent(X_DIM_GRID);
    int N_y = system.grid.extent(Y_DIM_GRID);
    int N_z = system.grid.extent(Z_DIM_GRID);
    // Create Functor
    const Initialize_grid_functor grid_initialize_functor(system);
    // Choose which initialization to use
    if (system.grid_initial_state == "testing") {
        if (system.print_output) {std::cout << "Initializing grid to a test configuration" << std::endl;}
        // Create Multi-Dim policy
        Kokkos::MDRangePolicy<grid_initialize_testing_tag, Kokkos::Rank<4>> policy(
                {0 + system.guardCells, 0 + system.guardCells, 0 + system.guardCells, 0},
                {N_z - system.guardCells, N_y - system.guardCells, N_x - system.guardCells, N_v});
        // Use Kokkos and functor to initialize grid
        Kokkos::parallel_for(policy, grid_initialize_functor);
    }
    else if (system.grid_initial_state == "advect_1d_hat"){
        if (system.print_output) {std::cout << "Initializing grid to 1D hat advection test" << std::endl;}
        // Create Multi-Dim policy
        // Note that I am filling the guard cells. This is needed in the 1D case.
        Kokkos::MDRangePolicy<grid_initialize_1d_advect_tag,Kokkos::Rank<4>> policy({0,0,0, 0},{N_z,N_y,N_x, N_v});
        // Use Kokkos and functor to initialize grid
        Kokkos::parallel_for(policy, grid_initialize_functor);
    }

    else{
        throw(std::invalid_argument){"Invalid grid_initial_state used."};
    }
}



#endif //HYDROSOLVER_INITIALIZE_GRID_H
