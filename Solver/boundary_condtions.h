//
// Created by carls502 on 1/20/2020.
//

#ifndef HYDROSOLVER_BOUNDARY_CONDTIONS_H
#define HYDROSOLVER_BOUNDARY_CONDTIONS_H



#include "System.h"


// Tag for different ways to apply boundary conditions
struct periodic_boundary_conditions{};


// ------------------- Functor for calculating dt given the system's state ------------- //
struct boundary_condition_functor{
    // TODO: have a different function fill left and right cells. Might help data locality
    // What dimension to apply the boundary conditions to
    // x = 0, y = 1, z = 2
    const int dim;
    const int num_guard_cells;
    const System system;

    explicit boundary_condition_functor(const System system, int num_guard_cells, int dim) : system(system), dim(dim), num_guard_cells(num_guard_cells){};

    // Kokkos function to apply boundary condition to.
    KOKKOS_INLINE_FUNCTION
    void operator() (const periodic_boundary_conditions&, const long unsigned int& g, const long unsigned int& j, const long unsigned int& i, const long unsigned int& v) const{
        /* Set periodic boundary conditions. The indices i and j correspond to either x, y, or z dim depending on which dim is set. It goes as follows:
         * x : i = y & j = z
         * y : i = z & j = x
         * z : i = x & j = y
         * See: IntroCompHydro 4.1
         */
        int master_R = system.grid.extent(dim) -  num_guard_cells*2 + g;
        int guard_L = g;
        int guard_R = system.grid.extent(dim) - num_guard_cells + g;
        int master_L = num_guard_cells + g;
        switch (dim){
            case X_DIM_GRID:
                system.grid(j, i, guard_L, v) = system.grid(j, i, master_R, v);
                system.grid(j, i, guard_R, v) = system.grid(j, i, master_L, v);
                break;
            case Y_DIM_GRID:
                system.grid(i, guard_R, j, v) = system.grid(i, master_L, j, v);
                system.grid(i, guard_L, j, v) = system.grid(i, master_R, j, v);
                break;
            case Z_DIM_GRID:
                system.grid(guard_L, j, i, v) = system.grid(master_R, j, i, v);
                system.grid(guard_R, j, i, v) = system.grid(master_L, j, i, v);
                break;
        }


    }

};

// ------------------ Callable Function for calculating time step ---------------------- //
void apply_boundary_conditions(const System& system){
    /*
     * Apply boundary conditions
     */
    // Get grid dimensions
    int N_v = system.grid.extent(VAR_DIM_GRID);
    int N_x = system.grid.extent(X_DIM_GRID);
    int N_y = system.grid.extent(Y_DIM_GRID);
    int N_z = system.grid.extent(Z_DIM_GRID);
    int guard_cells = system.guardCells;
    // apply boundary conditions to in each directions

    // X Dim
    if (system.simulation_dim >= 1){
        // Set x-dim
        auto bc = boundary_condition_functor(system, guard_cells, X_DIM_GRID);
        // Create range policy
        if (system.boundary_condition_x == "periodic") {
            Kokkos::MDRangePolicy<periodic_boundary_conditions, Kokkos::Rank<4>> policy({0, 0, 0, 0}, {guard_cells, N_z, N_y, N_v});
            Kokkos::parallel_for(policy, bc);
        }
        else {
            throw(std::invalid_argument){"Invalid boundary condition used for x dim."};
        }
    }if (system.simulation_dim >= 2){
        // Set y-dim
        auto bc = boundary_condition_functor(system, guard_cells, Y_DIM_GRID);
        // Create range policy
        if (system.boundary_condition_x == "periodic") {
            Kokkos::MDRangePolicy<periodic_boundary_conditions, Kokkos::Rank<4>> policy({0, 0, 0, 0}, {guard_cells, N_x, N_z, N_v});
            Kokkos::parallel_for(policy, bc);
        }
        else {
            throw(std::invalid_argument){"Invalid boundary condition used for y dim."};
        }
    }if (system.simulation_dim >= 3){
        // Set z-dim
        auto bc = boundary_condition_functor(system, guard_cells, Z_DIM_GRID);
        // Create range policy
        if (system.boundary_condition_x == "periodic") {
            Kokkos::MDRangePolicy<periodic_boundary_conditions, Kokkos::Rank<4>> policy({0, 0, 0, 0}, {guard_cells, N_y, N_x, N_v});
            Kokkos::parallel_for(policy, bc);
        }
        else {
            throw(std::invalid_argument){"Invalid boundary condition used for z dim."};
        }
    }

}



#endif //HYDROSOLVER_BOUNDARY_CONDTIONS_H
