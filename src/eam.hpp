#ifndef EAM_HPP_
#define EAM_HPP_

#include <memory>
#include <vector>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

namespace potential {
class EAM {
  public:
    EAM(double const rcut_, double const rhomax_fatal_, double const dr_interp_,
        double const drho_interp_, std::vector<double> const &rho_,
        std::vector<double> const &drhodr_, std::vector<double> const &d2rhodr2_,
        std::vector<double> const &F_, std::vector<double> const &dFdrho_,
        std::vector<double> const &d2Fdrho2_)
        : rcut(rcut_), rhomax_fatal(rhomax_fatal_),

          spline_rho(rho_.begin(), rho_.end(), 0., dr_interp_),
          spline_drhodr(drhodr_.begin(), drhodr_.end(), 0., dr_interp_),
          spline_d2rhodr2(d2rhodr2_.begin(), d2rhodr2_.end(), 0., dr_interp_),

          spline_F(F_.begin(), F_.end(), 0., drho_interp_),
          spline_dFdrho(dFdrho_.begin(), dFdrho_.end(), 0., drho_interp_),
          spline_d2Fdrho2(d2Fdrho2_.begin(), d2Fdrho2_.end(), 0., drho_interp_){};

    virtual ~EAM() = default;
    double const rcut;
    double const rhomax_fatal;

    boost::math::interpolators::cardinal_cubic_b_spline<double> const spline_rho;
    boost::math::interpolators::cardinal_cubic_b_spline<double> const spline_drhodr;
    boost::math::interpolators::cardinal_cubic_b_spline<double> const spline_d2rhodr2;

    boost::math::interpolators::cardinal_cubic_b_spline<double> const spline_F;
    boost::math::interpolators::cardinal_cubic_b_spline<double> const spline_dFdrho;
    boost::math::interpolators::cardinal_cubic_b_spline<double> const spline_d2Fdrho2;

    double rho(double const dr) {
        if (dr > rcut)
            return 0.;
        return spline_rho(dr);
    }

    double drhodr(double const dr) {
        if (dr > rcut)
            return 0.;
        return spline_drhodr(dr);
    }

    double d2rhodr2(double const dr) {
        if (dr > rcut)
            return 0.;
        return spline_d2rhodr2(dr);
    }

    double F(double const rho) {
        if (rho > rhomax_fatal) {
            printf("FATAL: EAM F rho >rhomax_fatal %g %g\n", rho, rhomax_fatal);
            exit(1);
        }
        return spline_F(rho);
    }

    double dFdrho(double const rho) {
        if (rho > rhomax_fatal) {
            printf("FATAL: EAM dFdrho rho >rhomax_fatal %g %g\n", rho, rhomax_fatal);
            exit(1);
        }
        return spline_dFdrho(rho);
    }

    double d2Fdrho2(double const rho) {
        if (rho > rhomax_fatal) {
            printf("FATAL: EAM d2Fdrho2 rho >rhomax_fatal %g %g\n", rho, rhomax_fatal);
            exit(1);
        }
        return spline_d2Fdrho2(rho);
    }
};
using EAMPtr = std::shared_ptr<EAM>;
} // namespace potential
#endif /* EAM_HPP_ */
