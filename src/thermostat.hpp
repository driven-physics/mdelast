#ifndef THERMOSTAT_HPP_
#define THERMOSTAT_HPP_

#include <vector>
#include <cmath>

//--- Nose-Hoover-Chain thermostat

struct Thermostat {
    double const temperature;    // Specified temperature
    std::vector<double> const q; // Thermal inertias
    std::vector<double> p_eta;   // Thermal momenta
    std::vector<double> eta;   // Thermal coordinates (needed only to calculate conserved quantity)
    unsigned int const nchain; // Number of Nose-Hoover chain variables
    double const Nf; // number of degrees of freedom for the calculation of thermostat energy;

    Thermostat(double const temperature_, std::vector<double> const &q_,
               std::vector<double> const &p_eta_, double const Nf_);

    void U4_propagator(double const dt,
                       double const e_couple, // the coupled energy
                       bool const outward);

    template <class T> void U3_propagator(double const dt, T *p_coupled, unsigned int size) {
        // Eq. is p -> p * exp(-dt * p_eta[0]/q[0]).
        // The factor is identical for all coupled DOFs, so compute once.
        double const scale = std::exp(-dt * p_eta[0] / q[0]);
        for (unsigned int i = 0; i < size; ++i) {
            p_coupled[i] *= scale;
        }

        for (unsigned int j = 0; j < nchain; ++j) {
            eta[j] += dt * p_eta[j] / q[j];
        }
    }

    double calculate_etot();
};

#endif /* THERMOSTAT_HPP_ */
