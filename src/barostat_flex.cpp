#include "barostat_flex.hpp"

BarostatFlex::BarostatFlex(double const pressure_, double const w_, Eigen::Matrix3d const pg_,
                           Thermostat const &thermostat_)
    : pressure(pressure_), w(w_), pg(pg_), thermostat(thermostat_) {
    assert(w > 0.);
}
