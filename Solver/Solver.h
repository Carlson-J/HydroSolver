//
// Created by carls502 on 1/13/2020.
//

#ifndef HYDROSOLVER_SOLVER_H
#define HYDROSOLVER_SOLVER_H

#include "../../kokkos/core/src/Kokkos_Core.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Initialize_grid.h"
#include "Calculate_dt.h"
#include "boundary_condtions.h"
#include "Solve_iteration.h"
#include "System.h"
#include <iomanip>
#include "Timer.h"
#include <json.hpp>

using json = nlohmann::json;

using std::string;
using std::cout;
using std::endl;
using std::string;

class Solver {

private:
    const System system;
    double t;
    double dt;
    int steps;
    int last_save_step;
    // Output file name. The timestep will be appended to it.
    string outfile_name;



public:
    Timer timer;


    explicit Solver(const System system, string outfile_name) : system(system), timer(){
        t = 0;
        dt = 0;
        steps = 0;
        last_save_step = -1;
        this->outfile_name = outfile_name;
    }


    // ----------- Simulation Driver -------------------------------------- //
    void Run_simulation(double runTime){
        // Runs simulation until runTime is reached. Starts at time 0
        // The overall methods is outlined in IntroCompHydro pg.52
        /*
         * Set initial conditions
         * Main evolution loop
         *  apply boundary conditions
         *  calculate time step
         *  compute interface states
         *  solve Riemann problem
         *  do conservative update
         *
         * Note that the last three are all done in Solve_iteration(), because their form will depend
         * on how we solve the time dependence, e.g. forward euler or Runge-kutta
         */

        // Save header
        if (system.intermediate_saving) {Save_header();}

        // Initialize grid
        Initialize_grid(system);

        // Reset time and steps
        t = 0;
        steps = 0;

        // Save initial solution
        if (system.intermediate_saving){
            if(system.print_output) {cout << "Saving Initial solution" << endl;}
            last_save_step = steps;
            Save_solution();
        }


        while (t < runTime){
//            // Print grid
//            for (int i = 0; i< system.grid.extent(Z_DIM_GRID); i++){
//                std::cout << "Z " << i << endl << toString_grid_values(2,i,0,true);
//            }
//            cout << endl << endl;

            // Print step information
            if (system.print_output){
                cout << toString_simulation_status();
            }
            // Apply Boundary Conditions
            apply_boundary_conditions(system);

            // Calculate time step
            dt = Calc_dt(runTime, t, system);

            // Solve iteration
            Solve_iteration(system);

            // Update time and steps
            Update_simulation_time();
            steps++;



            // Save solution step if enough steps have passed
            if (system.intermediate_saving && (system.save_frequency <= steps - last_save_step)){
                if(system.print_output) {cout << "Saving solution" << endl;}
                last_save_step = steps;
                Save_solution();
            }
        }
        if(system.print_output){
            cout << "Simulation finished after " << steps << " steps at time " << t << endl;
        }
    };


    // ----------- Functions called to run simulation --------------------- //
    void Update_simulation_time() { t += dt; }



    // ----------- Methods used to save system ---------------------------- //
    void Save_solution();
    void Save_header();


    // ---------- Diagnostic Methods -------------------------------------- //

    string toString_simulation_parameters();
    string toString_simulation_status();
    string toString_grid_values(int digits, int z, int var, bool print_guard_cells);
};





//------------------- Methods for getting printout of simulation Status -----------//
string Solver::toString_simulation_status(){
    /*
     * Returns status of simulation
     */
    std::ostringstream output;
    output << "Step: " << steps << " dt: " << dt << " Time: " << t << endl;
    return output.str();
}

string Solver::toString_simulation_parameters(){
    /*
     * Returns string of current simulation parameters
     */
    std::ostringstream output;
    output << "dt " << dt;
    return output.str();
}

string Solver::toString_grid_values(int digits, int z, int var, bool print_guard_cells){
    /*
     * Returns a string that shows the grid up to a number of digits. This will only print the first 2 dim.
     */
    int i0 = 0;
    int j0 = 0;
    int i_final = system.grid.extent(X_DIM_GRID);
    int j_final = system.grid.extent(Y_DIM_GRID);
    if (!print_guard_cells){
        i0 = system.guardCells;
        j0 = system.guardCells;
        i_final -= system.guardCells;
        j_final -= system.guardCells;
    }
    std::ostringstream output;
    output <<  std::scientific << std::setprecision(digits);
    for (int j = j0; j < j_final; j++) {
        for (int i = i0; i < i_final - 1; i++){
            output << std::setw(digits + 9) << system.grid(z,j,i,var) << ",";
        }
        output << std::setw(digits + 9) << system.grid(z, j, i_final - 1, var);
        output << endl;
    }
    return output.str();
}

void Solver::Save_solution(){
    // Save solution to file in binary format
    string stepsStr = std::to_string(steps);
    std::string new_string = outfile_name + "/data/" + std::string(6 - stepsStr.length(), '0') + stepsStr + ".sol";
    std::ofstream outfile( new_string,std::ofstream::out| std::ios_base::binary);
    /*
     * Save the grid array the same way it is laid out in memory.
     */
    if (!outfile.is_open()){
        throw(std::invalid_argument){"Data file could not be opened!!"};
    }
    for (int k =0 ;k < system.grid.extent(Z_DIM_GRID); k++) {
        for (int j = 0; j < system.grid.extent(Y_DIM_GRID); j++) {
            for (int i = 0; i < system.grid.extent(X_DIM_GRID); i++) {
                for (int v = 0; v < system.grid.extent(VAR_DIM_GRID); v++) {
                    outfile.write(reinterpret_cast<char*>(&system.grid(k,j,i,v)), sizeof(double));
//                    outfile << system.grid(v,i,j,k);
                }
            }
        }
    }

    outfile.close();

}

void Solver::Save_header(){
    // Save the header in a json file. This will give the structure of the arrays to be saved in binary
    std::ofstream header_file(outfile_name + "/header.json", std::ios_base::out);

    // Create json file to be saved as header.
    json header;
    header["Nx"] = system.grid.extent(X_DIM_GRID);
    header["Ny"] = system.grid.extent(Y_DIM_GRID);
    header["Nz"] = system.grid.extent(Z_DIM_GRID);
    header["Nv"] = system.grid.extent(VAR_DIM_GRID);
    header["GuardCells"] = system.guardCells;
    header["dx"] = system.dx;


    // Save file
    header_file << header.dump(3);


}

#endif //HYDROSOLVER_SOLVER_H
