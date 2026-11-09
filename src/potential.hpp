#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <boost/math/special_functions/pow.hpp>
#include "poly_table.hpp"
#include "units.hpp"

namespace potential {

using boost::math::pow;

class Potential {
  public:
    Potential(double const rcut_) : rcut(rcut_) {
        printf("Constructor of Potential rcut=%g [A]\n", rcut / Units::angstrom);
    }
    virtual ~Potential() {
        printf("Destructor of Potential\n");
    }
    virtual values calculate(double const dr) = 0;
    double const rcut;
};
using PotentialPtr = std::shared_ptr<Potential>;

class Morse : public Potential {
  public:
    Morse(double const rcut_, double const D_, double const alpha_, double const r0_)
        : Potential(rcut_), D(D_), alpha(alpha_), r0(r0_){};

    virtual ~Morse() = default;

    double const D;
    double const alpha;
    double const r0;

    values calculate(double const dr) override {
        if (dr > rcut) {
            return {0., 0., 0.};
        } else {
            double const exp1 = std::exp(-alpha * (dr - r0));
            double const exp2 = exp1 * exp1;

            return {D * (exp2 - 2. * exp1), -2. * D * alpha * (exp2 - exp1),
                    2. * alpha * alpha * D * (2 * exp2 - exp1)};
        }
    }
};

class LennardJones : public Potential {
  public:
    LennardJones(double const rcut_, double const epsilon_, double const sigma_)
        : Potential(rcut_), epsilon(epsilon_), sigma(sigma_),
          vcut(4. * epsilon_ * (pow<12>(sigma_ / rcut_) - pow<6>(sigma_ / rcut_))){};

    virtual ~LennardJones() = default;

    double const epsilon;
    double const sigma;
    double const vcut;

    values calculate(double const dr) override {
        if (dr > rcut) {
            return {0., 0., 0.};
        } else {
            double const v6 = pow<6>(sigma / dr);
            double const v12 = v6 * v6;
            return {4. * epsilon * (v12 - v6) - vcut, 24. * epsilon * (v6 - 2. * v12) / dr,
                    24. * epsilon * (26. * v12 - 7. * v6) / (dr * dr)};
        }
    }
};

class TablePotential : public Potential {
  public:
    TablePotential(double const rcut_, double const rmin_fatal_, double const dr_interp_,
                   std::vector<double> const &v_, std::vector<double> const &dvdr_,
                   std::vector<double> const &d2vdr2_,
                   std::string const &interpolation_ = "quintic")
        : Potential(rcut_), rmin_fatal(rmin_fatal_),
          table(interpolation_, dr_interp_, v_, dvdr_, d2vdr2_){};

    virtual ~TablePotential() = default;

    double const rmin_fatal;
    Table const table; // V, dV/dr, d2V/dr2 ("quintic": one piecewise quintic; "bspline": v1.1)

    values calculate(double const dr) override {

        if (dr < rmin_fatal) {
            printf("FATAL: PotentialTable dr <rmin_fatal %g %g\n", dr, rmin_fatal);
        }

        if (dr > rcut) {
            return {0., 0., 0.};

        } else {
            return table.all(dr);
        }
    }
};

} // namespace potential
#endif /* POTENTIAL_HPP_ */
