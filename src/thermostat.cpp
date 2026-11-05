#include "thermostat.hpp"
#include "units.hpp"
#include <cassert>
#include <stdio.h>
#include <boost/math/special_functions/pow.hpp>

#include "numerics.hpp"

using boost::math::pow;

Thermostat::Thermostat(double const temperature_, std::vector<double> const &q_,
                       std::vector<double> const &p_eta_, double const Nf_)
    : temperature(temperature_), q(q_), p_eta(p_eta_), eta(p_eta_.size(), 0.),
      nchain(p_eta_.size()), Nf(Nf_) {

    assert(temperature_ > 0.);
    assert(nchain == q.size());
    assert(nchain == p_eta.size());
    assert(nchain == eta.size());
    assert(nchain > 0);
    assert(Nf > 0.);
}

double Thermostat::calculate_etot() {
    //----- Thermostat part of Eq. 9 in: T.-Q. Yu et al. / Chemical Physics 370 (2010) 294305

    double etot = Nf * eta[0];
    for (unsigned int j = 1; j < nchain; ++j) {
        etot += eta[j];
    }
    etot *= Units::k_boltz * temperature;

    for (unsigned int j = 0; j < nchain; ++j) {
        etot += 0.5 * pow<2>(p_eta[j]) / q[j];
    }
    return etot;
}

//---- U4: thermostat propagator
void Thermostat::U4_propagator(double const dt,
                               double const e_couple, // the coupled energy
                               bool const outward) {

    int const dj = outward ? 1 : -1;
    int const j_first = outward ? 0 : nchain - 1;
    int const j_last = outward ? static_cast<int>(nchain) : -1;
    double const kT = Units::k_boltz * temperature;
    int const j_end = nchain - 1;

    for (int j = j_first; j != j_last; j += dj) {

        //----- Eq. 8 in: T.-Q. Yu et al. / Chemical Physics 370 (2010) 294305
        double const Gj = j != 0 ? pow<2>(p_eta[j - 1]) / q[j - 1] - kT : 2. * e_couple - Nf * kT;

        if (j != j_end) {

            //-- Allen-Tildesley eq. 3.74a
            double const x = dt * p_eta[j + 1] / q[j + 1];
            auto const res = numerics::expx_operator(x);
            p_eta[j] = p_eta[j] * res.expx + dt * Gj * res.op;
        } else {
            //-- Allen-Tildesley eq. 3.73
            p_eta[j] += dt * Gj;
        }
    }
}
