//
// Created by carls502 on 1/23/2020.
//

#ifndef HYDROSOLVER_TIMER_H
#define HYDROSOLVER_TIMER_H
#include <string>
#include <vector>
#include "../../kokkos/core/src/Kokkos_Core.hpp"
using std::string;
using std::vector;

Kokkos::Timer kokkos_timer;
struct Event{
//    double duration = 0;
    double start_time = 0;
    double end_time = -1;
    string label;

    Event(string event_name){
        label = std::move(event_name);
        start_time = kokkos_timer.seconds();
    }
    double Get_time(){
        if (end_time == -1){
            throw std::runtime_error("Time was not stopped before Get_time was called");
        }
        return end_time - start_time;
    }
};

class Timer{
public:
    void Start_timer(std::string event_name){
        // Creates a new timer or resets one if the label already exists.
        int id = Find_event(event_name);
        if (id == -1){
            events.emplace_back(Event(std::move(event_name)));
        }
        else {
            events.at(id).start_time = kokkos_timer.seconds();
        }

    }

    void Stop_time(string event_name){
        // Resets timer and saves duration from last reset or creation.
        double time = kokkos_timer.seconds();
        int id = Find_event(event_name);
        if (id < 0) {
            throw std::invalid_argument("Event has not been created, but Stop_time has been called");
        }
        events.at(id).end_time = time;
    }

    double Get_time(string event_name){
        // Get the duration of a timer that has been stopped.
        int id = Find_event(event_name);
        if (id < 0) {
            throw std::invalid_argument("Event has not been created, but Get_time has been called");
        }
        return events.at(id).Get_time();
    }

private:
    vector<Event> events;
    int Find_event(string event_name){
        // Find a timer with a specific label. Return the index to the timer or -1 if it does not exist.
        for (int i = 0; i < events.size(); i++){
            if (events.at(i).label == event_name){
                return i;
            }
        }
        return -1;
    }

};

#endif //HYDROSOLVER_TIMER_H
