// Unit tests of potential::QuinticTable (src/poly_table.hpp).
//
// Build and run (repository warning flags plus the address and undefined-behaviour
// sanitizers, which turn any out-of-range table access into a hard error):
//   c++ -std=c++17 -O1 -g -Wall -Wextra -Werror -pedantic-errors -Wshadow -Wdouble-promotion
//     -fsanitize=address,undefined -fno-sanitize-recover=all -I src -I <boost include>
//     tools/test_poly_table.cpp -o build/test_poly_table && build/test_poly_table
//
// Tests
//   1 polynomial exactness   every polynomial of degree <= 5 is reproduced (f, f', f'')
//   2 convergence order      smooth f: error ~ h^6, h^5, h^4 for f, f', f''
//   3 node values            f, f', f'' at the nodes equal the input arrays
//   4 continuity             left/right limits at every interior node, derivatives 0..5
//                            (0..2 must agree to roundoff for ANY input, 3..5 need not)
//   5 C2 for arbitrary data  continuity 0..2 with random (inconsistent) f, f', f''
//   6 index edge cases       x = 0, -0, last node, nextafter around nodes, extrapolation
//                            below 0 and beyond the end, huge |x|, +-inf, NaN
//   7 minimal table          two nodes (one interval)
// Exit status 0 when every check passes.
#include "poly_table.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <limits>
#include <random>
#include <string>
#include <vector>

using potential::QuinticTable;
using ld = long double;

static int n_fail = 0, n_pass = 0;

static void check(bool ok, std::string const &what, double val, double tol) {
    (ok ? n_pass : n_fail)++;
    printf("  %s  %-58s %10.3e (tol %.1e)\n", ok ? "PASS" : "FAIL", what.c_str(), val, tol);
}

struct Fun {
    std::function<double(double)> f, f1, f2;
};

static QuinticTable make(Fun const &F, double h, int n) {
    std::vector<double> a(n + 1), b(n + 1), c(n + 1);
    for (int i = 0; i <= n; ++i) {
        double x = i * h;
        a[i] = F.f(x), b[i] = F.f1(x), c[i] = F.f2(x);
    }
    return QuinticTable(h, a, b, c);
}

// derivative k (0..5) of the polynomial of interval i at local coordinate t, in long double
static ld poly_deriv(QuinticTable const &T, int i, int k, ld t) {
    auto const &c = T.coef[static_cast<std::size_t>(i)];
    ld s = 0;
    for (int m = 5; m >= k; --m) {
        ld fac = 1;
        for (int j = 0; j < k; ++j)
            fac *= (m - j);
        s = s * t + fac * c[m];
    }
    return s / std::pow(static_cast<ld>(T.h), k);
}

// max relative error of f, f', f'' of table T against F over [lo, hi]
static void max_err(QuinticTable const &T, Fun const &F, double lo, double hi, double out[3]) {
    double sc[3] = {0, 0, 0}, e[3] = {0, 0, 0};
    for (int j = 0; j <= 20000; ++j) {
        double x = lo + (hi - lo) * j / 20000.;
        double r[3] = {F.f(x), F.f1(x), F.f2(x)};
        double v[3] = {T.value(x), T.d1(x), T.d2(x)};
        auto a = T.all(x);
        double w[3] = {a.V, a.dVdr, a.d2Vdr2};
        for (int k = 0; k < 3; ++k) {
            sc[k] = std::fmax(sc[k], std::fabs(r[k]));
            e[k] = std::fmax(e[k], std::fabs(v[k] - r[k]));
            if (v[k] != w[k]) // all() must be bit-identical to value/d1/d2
                e[k] = INFINITY;
        }
    }
    for (int k = 0; k < 3; ++k)
        out[k] = e[k] / std::fmax(sc[k], 1e-300);
}

int main() {
    char buf[256];

    // Input rounding limits the attainable accuracy between nodes: a perturbation df of a
    // node value changes f' by up to 1.875 df/h and f'' by up to 5.774 df/h^2 (maxima of the
    // first and second derivatives of the quintic Hermite value basis 1 - 10t^3 + 15t^4 - 6t^5).
    // With df = 4 ulp of max|f| the bounds below are the conditioning limit of ANY interpolant
    // through rounded node values, not a property of this implementation.
    printf("1. polynomials of degree 0..5 are reproduced exactly\n");
    for (double const h : {0.1, 0.0034}) {
        int const n = static_cast<int>(std::lround(6.8 / h));
        for (int deg = 0; deg <= 5; ++deg) {
            // p(x) = sum_m (x - 3)^m / (m + 1), centred so that the terms do not cancel trivially
            Fun P{[deg](double x) { double s = 0; for (int m = deg; m >= 0; --m) s = s * (x - 3) + 1. / (m + 1); return s; },
                  [deg](double x) { double s = 0; for (int m = deg; m >= 1; --m) s = s * (x - 3) + m / (m + 1.); return s; },
                  [deg](double x) { double s = 0; for (int m = deg; m >= 2; --m) s = s * (x - 3) + m * (m - 1) / (m + 1.); return s; }};
            auto T = make(P, h, n);
            double e[3];
            max_err(T, P, 0., n * h, e);
            double fmax = 0, sc[3] = {0, 0, 0};
            for (int j = 0; j <= 2000; ++j) {
                double x = n * h * j / 2000.;
                fmax = std::fmax(fmax, std::fabs(P.f(x)));
                sc[1] = std::fmax(sc[1], std::fabs(P.f1(x)));
                sc[2] = std::fmax(sc[2], std::fabs(P.f2(x)));
            }
            double const df = 4 * 2.22e-16 * fmax;
            double const tol[3] = {1e-14, std::fmax(1e-13, 2 * 1.875 * df / h / sc[1]),
                                   std::fmax(1e-13, 2 * 5.774 * df / (h * h) / std::fmax(sc[2], 1e-300))};
            for (int k = 0; k <= std::min(deg, 2); ++k) {
                snprintf(buf, sizeof buf, "h %.4g, degree %d, derivative %d, max rel err", h, deg, k);
                check(e[k] <= tol[k], buf, e[k], tol[k]);
            }
        }
    }

    printf("2. convergence order on f = sin(3x) exp(-x/2) + 1/(1+x^2) over [0, 4]\n");
    Fun S{[](double x) { return std::sin(3 * x) * std::exp(-x / 2) + 1 / (1 + x * x); },
          [](double x) { return std::exp(-x / 2) * (3 * std::cos(3 * x) - 0.5 * std::sin(3 * x)) - 2 * x / ((1 + x * x) * (1 + x * x)); },
          [](double x) { return std::exp(-x / 2) * (-8.75 * std::sin(3 * x) - 3 * std::cos(3 * x)) + (6 * x * x - 2) / std::pow(1 + x * x, 3); }};
    double prev[3] = {0, 0, 0};
    int const expect[3] = {6, 5, 4};
    for (int n : {20, 40, 80, 160}) {
        auto T = make(S, 4. / n, n);
        double e[3];
        max_err(T, S, 0., 4., e);
        printf("   n = %4d  err f %.3e  f' %.3e  f'' %.3e\n", n, e[0], e[1], e[2]);
        if (n > 20)
            for (int k = 0; k < 3; ++k) {
                double order = std::log2(prev[k] / e[k]);
                snprintf(buf, sizeof buf, "n %d: observed order of derivative %d (expect %d)", n, k, expect[k]);
                check(std::fabs(order - expect[k]) < 0.35, buf, order, 0.35);
            }
        for (int k = 0; k < 3; ++k) prev[k] = e[k];
    }

    printf("3. node values equal the input arrays (n = 2000)\n");
    {
        int const n = 2000;
        double const h = 0.0034;
        auto T = make(S, h, n);
        double e[3] = {0, 0, 0};
        for (int i = 0; i <= n; ++i) {
            double x = i * h;
            e[0] = std::fmax(e[0], std::fabs(T.value(x) - S.f(x)) / 2.);
            e[1] = std::fmax(e[1], std::fabs(T.d1(x) - S.f1(x)) / 4.);
            e[2] = std::fmax(e[2], std::fabs(T.d2(x) - S.f2(x)) / 12.);
        }
        check(e[0] < 1e-14, "f   at nodes (rel)", e[0], 1e-14);
        check(e[1] < 1e-12, "f'  at nodes (rel)", e[1], 1e-12);
        check(e[2] < 1e-10, "f'' at nodes (rel)", e[2], 1e-10);
    }

    printf("4/5. continuity at interior nodes: |left - right| / max|f^(k)|\n");
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> U(-1, 1);
    for (int pass = 0; pass < 2; ++pass) {
        int const n = 400;
        double const h = 4. / n;
        QuinticTable T = [&] {
            if (pass == 0) return make(S, h, n);
            std::vector<double> a(n + 1), b(n + 1), c(n + 1);
            for (int i = 0; i <= n; ++i) a[i] = U(rng), b[i] = U(rng) / h, c[i] = U(rng) / (h * h);
            return QuinticTable(h, a, b, c);
        }();
        for (int k = 0; k <= 5; ++k) {
            ld jmax = 0, sc = 0;
            for (int i = 1; i < n; ++i) {
                ld L = poly_deriv(T, i - 1, k, 1), R = poly_deriv(T, i, k, 0);
                jmax = std::fmax(jmax, std::fabs(L - R));
                sc = std::fmax(sc, std::fmax(std::fabs(L), std::fabs(R)));
            }
            double rel = static_cast<double>(jmax / sc);
            snprintf(buf, sizeof buf, "%s data, derivative %d jump", pass ? "random" : "smooth", k);
            if (k <= 2)
                check(rel < (k == 2 ? 1e-12 : 1e-13), buf, rel, k == 2 ? 1e-12 : 1e-13);
            else
                printf("  info  %-58s %10.3e (C%d not imposed)\n", buf, rel, k);
        }
        // the same through the public interface, with x one ulp on either side of each node
        double e[3] = {0, 0, 0}, s[3] = {0, 0, 0};
        for (int i = 1; i < n; ++i) {
            double x = i * h, xl = std::nextafter(x, -1.), xr = std::nextafter(x, 10.);
            double l[3] = {T.value(xl), T.d1(xl), T.d2(xl)}, r[3] = {T.value(xr), T.d1(xr), T.d2(xr)};
            for (int k = 0; k < 3; ++k) {
                e[k] = std::fmax(e[k], std::fabs(l[k] - r[k]));
                s[k] = std::fmax(s[k], std::fabs(l[k]));
            }
        }
        for (int k = 0; k < 3; ++k) {
            snprintf(buf, sizeof buf, "%s data, x_i -+ 1 ulp, derivative %d (public API)", pass ? "random" : "smooth", k);
            check(e[k] / s[k] < 1e-11, buf, e[k] / s[k], 1e-11);
        }
    }

    printf("6. index edge cases (n = 2000, h = 0.0034, x_max = 6.8)\n");
    {
        int const n = 2000;
        double const h = 0.0034, xmax = n * h;
        auto T = make(S, h, n);
        auto finite3 = [&](double x) { return std::isfinite(T.value(x)) && std::isfinite(T.d1(x)) && std::isfinite(T.d2(x)); };
        // exact ends
        check(std::fabs(T.value(0.) - S.f(0.)) < 1e-15, "x = 0 gives f(0)", std::fabs(T.value(0.) - S.f(0.)), 1e-15);
        check(T.value(-0.) == T.value(0.), "x = -0.0 same as +0.0", 0, 0);
        check(std::fabs(T.value(xmax) - S.f(xmax)) < 1e-14, "x = x_max (last node) gives f(x_max)", std::fabs(T.value(xmax) - S.f(xmax)), 1e-14);
        double xm1 = std::nextafter(xmax, 0.), xp1 = std::nextafter(xmax, 1e9);
        check(std::fabs(T.value(xm1) - T.value(xp1)) < 1e-14, "x_max -+ 1 ulp continuous", std::fabs(T.value(xm1) - T.value(xp1)), 1e-14);
        // extrapolation = the end polynomial continued (no jump at the ends)
        double lo = std::nextafter(0., -1.);
        check(std::fabs(T.d2(lo) - T.d2(0.)) < 1e-10, "x = -1 ulp: f'' continuous with x = 0", std::fabs(T.d2(lo) - T.d2(0.)), 1e-10);
        double xe = xmax + 0.5 * h;
        ld ref = poly_deriv(T, n - 1, 0, 1.5L);
        check(std::fabs(T.value(xe) - static_cast<double>(ref)) < 1e-14, "x = x_max + h/2 uses last polynomial", std::fabs(T.value(xe) - static_cast<double>(ref)), 1e-14);
        // sweep every node, one ulp below/at/above: all finite, interval index in range
        bool ok = true;
        for (int i = 0; i <= n; ++i)
            for (double x : {std::nextafter(i * h, -1.), i * h, std::nextafter(i * h, 1e9)})
                ok = ok && finite3(x);
        check(ok, "all nodes -1 ulp / exact / +1 ulp finite", 0, 0);
        // pathological arguments: must not read outside the table (ASan/UBSan abort otherwise)
        for (double x : {-1e300, -1e10, 1e10, 3e9 * h, 1e300}) {
            double v = T.value(x) + T.d1(x) + T.d2(x); // reaching the next line = no invalid read
            snprintf(buf, sizeof buf, "x = %g: no out-of-range access (sanitizers)", x);
            check(true, buf, v, 0);
        }
        for (double x : {INFINITY, -INFINITY, NAN}) {
            double v = T.value(x);
            snprintf(buf, sizeof buf, "x = %g: no out-of-range access, result non-finite", x);
            check(!std::isfinite(v), buf, v, 0);
        }
    }

    printf("7. two-node table (one interval) reproduces a quintic\n");
    {
        Fun P{[](double x) { return 1 + x * (2 + x * (-1 + x * (0.5 + x * (0.25 - 0.125 * x)))); },
              [](double x) { return 2 + x * (-2 + x * (1.5 + x * (1 - 0.625 * x))); },
              [](double x) { return -2 + x * (3 + x * (3 - 2.5 * x)); }};
        auto T = make(P, 1.0, 1);
        double e[3];
        max_err(T, P, 0., 1., e);
        check(std::fmax(e[0], std::fmax(e[1], e[2])) < 1e-14, "f, f', f'' max rel err", std::fmax(e[0], std::fmax(e[1], e[2])), 1e-14);
    }

    printf("\n=== %d passed, %d failed\n", n_pass, n_fail);
    return n_fail ? 1 : 0;
}
