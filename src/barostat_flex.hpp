#ifndef BAROSTAT_FLEX_HPP_
#define BAROSTAT_FLEX_HPP_

#include <Eigen/Dense>

#include "thermostat.hpp"

struct BarostatFlex {
    double const pressure; // Specified pressure
    double const w;        // barostat inertia
    Eigen::Matrix3d pg;
    Thermostat thermostat; //-- the NHC thermostat which is coupled to this barostat

    BarostatFlex(double const pressure_, double const w_, Eigen::Matrix3d const pg_,
                 Thermostat const &thermostat_);

    double calculate_etot();
};

#endif /* BAROSTAT_FLEX_HPP_ */
