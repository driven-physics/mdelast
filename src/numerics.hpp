#ifndef NUMERICS_HPP_
#define NUMERICS_HPP_
#include <array>

namespace numerics {

struct expx_result {
    double const expx;
    double const op;
};
static expx_result expx_operator(double const x) {
    // 5th order taylor expansion of (1-exp(-x))/x
    static double const c1 = -1.0 / 2.0, c2 = 1.0 / 6.0, c3 = -1.0 / 24.0, c4 = 1. / 120.,
                        c5 = -1. / 720.;
    double const expx = std::exp(-x);
    return {expx, std::abs(x) > 1e-3 ? (1.0 - expx) / x
                                     : 1.0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * c5))))};
}

// Seven-stage, sixth-order Suzuki-Yoshida coefficients
//	taken from Tuckerman book eq. 4.11.12
static std::array<double, 7> constexpr wsy7 = {
    0.784513610477560, 0.235573213359357,
    -1.17767998417887, 1. - 2. * (0.784513610477560 + 0.235573213359357 - 1.17767998417887),
    -1.17767998417887, 0.235573213359357,
    0.784513610477560,
};

} // namespace numerics

#endif /* NUMERICS_HPP_ */