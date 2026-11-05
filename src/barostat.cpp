#include "barostat.hpp"
#include <boost/math/special_functions/pow.hpp>

using boost::math::pow;

Barostat::Barostat(double const pressure_, double const w_, double const p_epsilon_,
                   Thermostat const &thermostat_)
    : pressure(pressure_), w(w_), epsilon(0.), p_epsilon(p_epsilon_), thermostat(thermostat_) {
    assert(w > 0.);
}

double Barostat::calculate_etot() {

    //----- Barostat part of part of Eq. 9 in: T.-Q. Yu et al. / Chemical Physics 370 (2010) 294305
    return 0.5 * pow<2>(p_epsilon) / w + thermostat.calculate_etot();
}
