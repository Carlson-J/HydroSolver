//
// Created by carls502 on 1/29/2020.
//

#ifndef HYDROSOLVER_ITERATION_SPACE_H
#define HYDROSOLVER_ITERATION_SPACE_H

#include "System.h"

struct Iteration_space{
    // Starting and finishing iteration indices for the gird.
    int N_v_0 = 0;
    int N_v_f = 0;
    int N_x_0 = 0;
    int N_x_f = 0;
    int N_y_0 = 0;
    int N_y_f = 0;
    int N_z_0 = 0;
    int N_z_f = 0;
    // Starting and ending iteration indices for the grid_flux
    int N_flux_v_0 = 0;
    int N_flux_v_f = 0;
    int N_flux_x_0 = 0;
    int N_flux_x_f = 0;
    int N_flux_y_0 = 0;
    int N_flux_y_f = 0;
    int N_flux_z_0 = 0;
    int N_flux_z_f = 0;

    explicit Iteration_space(const System system){
        initialize_grid_iteration_space(system);
        initialize_flux_grid_iteration_space(system);

    }

private:
    void initialize_grid_iteration_space(const System &system) {// initialize to the iteration space for a given dimension and number of guard cells.
        N_v_0 = 0;
        N_v_f = system.grid.extent(VAR_DIM_GRID);
        N_x_0 = 0;
        N_x_f = system.grid.extent(X_DIM_GRID);
        N_y_0 = 0;
        N_y_f = system.grid.extent(Y_DIM_GRID);
        N_z_0 = 0;
        N_z_f = system.grid.extent(Z_DIM_GRID);
        switch (system.simulation_dim) { // Note that there are no breaks, so the all cases below the one that triggers run.
            case 3:
                N_z_0 += system.guardCells;
                N_z_f -= system.guardCells;
            case 2:
                N_y_0 += system.guardCells;
                N_y_f -= system.guardCells;
            case 1:
                N_x_0 += system.guardCells;
                N_x_f -= system.guardCells;
        }
    }

    void initialize_flux_grid_iteration_space(const System &system) {// initialize to the iteration space for a given dimension and number of guard cells.
        N_flux_v_0 = 0;
        N_flux_v_f = system.grid_flux.extent(VAR_DIM_FLUX);
        N_flux_x_0 = 0;
        N_flux_x_f = system.grid_flux.extent(X_DIM_FLUX);
        N_flux_y_0 = 0;
        N_flux_y_f = system.grid_flux.extent(Y_DIM_FLUX);
        N_flux_z_0 = 0;
        N_flux_z_f = system.grid_flux.extent(Z_DIM_FLUX);
        switch (system.simulation_dim) { // Note that there are no breaks, so the all cases below the one that triggers run.
            case 3:
                N_flux_z_0 += system.guardCells - 1;
                N_flux_z_f -= system.guardCells - 1;
            case 2:
                N_flux_y_0 += system.guardCells - 1;
                N_flux_y_f -= system.guardCells - 1;
            case 1:
                N_flux_x_0 += system.guardCells - 1;
                N_flux_x_f -= system.guardCells - 1;
        }
    }
};
#endif //HYDROSOLVER_ITERATION_SPACE_H
