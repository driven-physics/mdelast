#ifndef POLY_TABLE_HPP_
#define POLY_TABLE_HPP_

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

namespace potential {

struct values {
    double const V;
    double const dVdr;
    double const d2Vdr2;
};

// Piecewise quintic table on a uniform grid x_i = i*h, i = 0..n.
//
// On every interval [x_i, x_{i+1}] the function is the quintic Hermite polynomial fixed by
// f, f' and f'' at the two nodes. The value, first and second derivatives are evaluated from
// the same polynomial, so f' and f'' are the exact derivatives of the interpolant, and f,
// f', f'' are continuous across the nodes for any input (C2). When the node data come from
// a quintic spline whose knots are grid nodes (scipy UnivariateSpline(k=5, s=0), as in
// get_potential.py), the table reproduces that spline exactly and is C4.
//
// Each interval stores 6 coefficients in t = (x - x_i)/h in [0, 1], so one lookup is one
// index computation, one contiguous 48-byte load, and a Horner evaluation (cf. the cubic
// per-interval coefficients of LAMMPS pair_eam, which are C1).
class QuinticTable {
  public:
    QuinticTable(double const h_, std::vector<double> const &f, std::vector<double> const &f1,
                 std::vector<double> const &f2)
        : h(h_), inv_h(1. / h_), inv_h2(1. / (h_ * h_)),
          n(f.size() > 0 ? static_cast<int>(f.size()) - 1 : 0), coef(n > 0 ? n : 0) {
        if (f.size() < 2 or f1.size() != f.size() or f2.size() != f.size() or not(h > 0.)) {
            printf("FATAL: QuinticTable needs >= 2 nodes, equal array sizes, and h > 0 "
                   "(sizes %zu %zu %zu, h %g)\n",
                   f.size(), f1.size(), f2.size(), h);
            exit(1);
        }
        // a3..a5 are O(h^3 f''') but result from cancellation of O(h f') terms; form them in
        // extended precision so that the stored coefficients, and therefore f'', are accurate.
        using ld = long double;
        ld const hl = static_cast<ld>(h);
        for (int i = 0; i < n; ++i) {
            std::size_t const a = static_cast<std::size_t>(i), b = a + 1;
            ld const dy = static_cast<ld>(f[b]) - static_cast<ld>(f[a]);
            ld const D0 = hl * static_cast<ld>(f1[a]), D1 = hl * static_cast<ld>(f1[b]);
            ld const S0 = hl * hl * static_cast<ld>(f2[a]), S1 = hl * hl * static_cast<ld>(f2[b]);
            coef[a] = {f[a],
                       static_cast<double>(D0),
                       static_cast<double>(0.5L * S0),
                       static_cast<double>(10.L * dy - 6.L * D0 - 4.L * D1 -
                                           0.5L * (3.L * S0 - S1)),
                       static_cast<double>(-15.L * dy + 8.L * D0 + 7.L * D1 +
                                           0.5L * (3.L * S0 - 2.L * S1)),
                       static_cast<double>(6.L * dy - 3.L * (D0 + D1) - 0.5L * (S0 - S1))};
        }
    }

    double const h;
    double const inv_h;
    double const inv_h2;
    int const n; // number of intervals
    std::vector<std::array<double, 6>> coef;

    double value(double const x) const {
        double t;
        auto const &c = locate(x, t);
        return c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * (c[4] + t * c[5]))));
    }

    double d1(double const x) const {
        double t;
        auto const &c = locate(x, t);
        return (c[1] + t * (2. * c[2] + t * (3. * c[3] + t * (4. * c[4] + t * 5. * c[5])))) *
               inv_h;
    }

    double d2(double const x) const {
        double t;
        auto const &c = locate(x, t);
        return (2. * c[2] + t * (6. * c[3] + t * (12. * c[4] + t * 20. * c[5]))) * inv_h2;
    }

    values all(double const x) const {
        double t;
        auto const &c = locate(x, t);
        return {c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * (c[4] + t * c[5])))),
                (c[1] + t * (2. * c[2] + t * (3. * c[3] + t * (4. * c[4] + t * 5. * c[5])))) *
                    inv_h,
                (2. * c[2] + t * (6. * c[3] + t * (12. * c[4] + t * 20. * c[5]))) * inv_h2};
    }

  private:
    // Interval containing x; t is the local coordinate. Points outside [0, n*h] use the end
    // polynomials (extrapolation), like the B-splines this replaces.
    // The range test is done in floating point before the conversion to int: converting a
    // value outside the int range (or NaN) is undefined behaviour and, on x86, yields INT_MIN,
    // i.e. a read far outside the table. A NaN argument takes i = 0 and returns NaN.
    std::array<double, 6> const &locate(double const x, double &t) const {
        double const u = x * inv_h;
        int i = 0;
        if (u >= static_cast<double>(n))
            i = n - 1;
        else if (u > 0.)
            i = static_cast<int>(u);
        t = u - static_cast<double>(i);
        return coef[static_cast<std::size_t>(i)];
    }
};

// The tabulation used up to version 1.1: three independent Boost cubic B-splines through the
// node values of f, f' and f'' (f' and f'' are not the derivatives of the f interpolant).
// Kept, unchanged, so that results of earlier versions can be reproduced.
class BSplineTable {
  public:
    BSplineTable(double const h, std::vector<double> const &f, std::vector<double> const &f1,
                 std::vector<double> const &f2)
        : s0(f.begin(), f.end(), 0., h), s1(f1.begin(), f1.end(), 0., h),
          s2(f2.begin(), f2.end(), 0., h) {}

    boost::math::interpolators::cardinal_cubic_b_spline<double> const s0;
    boost::math::interpolators::cardinal_cubic_b_spline<double> const s1;
    boost::math::interpolators::cardinal_cubic_b_spline<double> const s2;

    double value(double const x) const { return s0(x); }
    double d1(double const x) const { return s1(x); }
    double d2(double const x) const { return s2(x); }
    values all(double const x) const { return {s0(x), s1(x), s2(x)}; }
};

// Interpolation scheme of a tabulated function, chosen when the potential is loaded:
//   "quintic" (default)  QuinticTable  - piecewise quintic Hermite, consistent derivatives
//   "bspline"            BSplineTable  - the cubic B-splines of versions <= 1.1
class Table {
  public:
    enum class Kind { quintic, bspline };

    static Kind parse(std::string const &name) {
        if (name == "quintic")
            return Kind::quintic;
        if (name == "bspline")
            return Kind::bspline;
        printf("FATAL: unknown table interpolation '%s' (use 'quintic' or 'bspline')\n",
               name.c_str());
        exit(1);
    }

    Table(std::string const &kind_, double const h, std::vector<double> const &f,
          std::vector<double> const &f1, std::vector<double> const &f2)
        : kind(parse(kind_)), name(kind_), quintic(h, f, f1, f2),
          bspline(kind == Kind::bspline ? std::make_unique<BSplineTable>(h, f, f1, f2) : nullptr) {}

    Kind const kind;
    std::string const name;
    QuinticTable const quintic;
    std::unique_ptr<BSplineTable const> const bspline;

    double value(double const x) const {
        return kind == Kind::quintic ? quintic.value(x) : bspline->value(x);
    }
    double d1(double const x) const {
        return kind == Kind::quintic ? quintic.d1(x) : bspline->d1(x);
    }
    double d2(double const x) const {
        return kind == Kind::quintic ? quintic.d2(x) : bspline->d2(x);
    }
    values all(double const x) const {
        return kind == Kind::quintic ? quintic.all(x) : bspline->all(x);
    }
};

} // namespace potential
#endif /* POLY_TABLE_HPP_ */
