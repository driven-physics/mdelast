#ifndef EAM_HPP_
#define EAM_HPP_

#include <memory>
#include <string>
#include <vector>
#include "poly_table.hpp"

namespace potential {
class EAM {
  public:
    EAM(double const rcut_, double const rhomax_fatal_, double const dr_interp_,
        double const drho_interp_, std::vector<double> const &rho_,
        std::vector<double> const &drhodr_, std::vector<double> const &d2rhodr2_,
        std::vector<double> const &F_, std::vector<double> const &dFdrho_,
        std::vector<double> const &d2Fdrho2_, std::string const &interpolation_ = "quintic")
        : rcut(rcut_), rhomax_fatal(rhomax_fatal_),
          table_rho(interpolation_, dr_interp_, rho_, drhodr_, d2rhodr2_),
          table_F(interpolation_, drho_interp_, F_, dFdrho_, d2Fdrho2_){};

    virtual ~EAM() = default;
    double const rcut;
    double const rhomax_fatal;

    // rho(r) and F(rho) with first and second derivatives ("quintic": one piecewise quintic
    // per function; "bspline": three cubic B-splines per function, as in v1.1)
    Table const table_rho;
    Table const table_F;

    double rho(double const dr) {
        if (dr > rcut)
            return 0.;
        return table_rho.value(dr);
    }

    double drhodr(double const dr) {
        if (dr > rcut)
            return 0.;
        return table_rho.d1(dr);
    }

    double d2rhodr2(double const dr) {
        if (dr > rcut)
            return 0.;
        return table_rho.d2(dr);
    }

    double F(double const rho) {
        if (rho > rhomax_fatal) {
            printf("FATAL: EAM F rho >rhomax_fatal %g %g\n", rho, rhomax_fatal);
            exit(1);
        }
        return table_F.value(rho);
    }

    double dFdrho(double const rho) {
        if (rho > rhomax_fatal) {
            printf("FATAL: EAM dFdrho rho >rhomax_fatal %g %g\n", rho, rhomax_fatal);
            exit(1);
        }
        return table_F.d1(rho);
    }

    double d2Fdrho2(double const rho) {
        if (rho > rhomax_fatal) {
            printf("FATAL: EAM d2Fdrho2 rho >rhomax_fatal %g %g\n", rho, rhomax_fatal);
            exit(1);
        }
        return table_F.d2(rho);
    }
};
using EAMPtr = std::shared_ptr<EAM>;
} // namespace potential
#endif /* EAM_HPP_ */
