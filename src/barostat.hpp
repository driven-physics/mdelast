#ifndef BAROSTAT_HPP_
#define BAROSTAT_HPP_

#include "thermostat.hpp"

struct Barostat {
    double const pressure; // Specified pressure
    double const w;        // barostat inertia
    double epsilon;
    double p_epsilon;
    Thermostat thermostat; //-- the NHC thermostat which is coupled to this barostat

    Barostat(double const pressure_, double const w_, double const p_epsilon_,
             Thermostat const &thermostat_);

    double calculate_etot();
};

#endif /* BAROSTAT_HPP_ */
