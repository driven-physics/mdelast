#ifndef UNITS_HPP_
#define UNITS_HPP_

namespace Units {
class BindingClass {};
static double constexpr k_boltz = 1.38064852e-16;
static double constexpr Navogadro = 6.02214086e23; // mol-1

static double constexpr ev = 1.6021773e-12;
static double constexpr joule = 1e7;

static double constexpr angstrom = 1e-8;
static double constexpr nm = 1e-7;

static double constexpr bar = 1e6;
static double constexpr kbar = 1e3 * bar;
static double constexpr atm = 1.01325 * bar;
static double constexpr Pa = 1e-5 * bar;
static double constexpr GPa = 1e9 * Pa;

static double constexpr mus = 1e-6;
static double constexpr ns = 1e-9;
static double constexpr ps = 1e-12;
static double constexpr fs = 1e-15;
} // namespace Units

#endif /* UNITS_HPP_ */
