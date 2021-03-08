#ifndef __DERIVATIVES_H_INCLUDED__
#define __DERIVATIVES_H_INCLUDED__

#include <cassert>

namespace diff {
typedef double (*Func)(double);

inline double GetDerivativeFirstOrderDirect(Func f, double at, double step) {
    assert(step != 0);

    return (f(at + step) - f(at)) / step;
}

inline double GetDerivativeFirstOrderInverse(Func f, double at, double step) {
    assert(step != 0);

    return (f(at) - f(at - step)) / step;
}

inline double GetDerivativeFirstOrderSimmetrical(Func f, double at, double step) {
    assert(step != 0);

    return (f(at + step) - f(at - step)) / (2 * step);
}

inline double GetDerivativeSecondOrder(Func f, double at, double step) {
    assert(step != 0);

    return (4.0 / 3.0) * (f(at + step) - f(at - step)) / (2 * step) -
           (1.0 / 3.0) * (f(at + 2 * step) - f(at - 2 * step)) / (4 * step);
}

inline double GetDerivativeThirdOrder(Func f, double at, double step) {
    assert(step != 0);

    return (3.0 / 2.0) * (f(at + step) - f(at - step)) / (2 * step) -
           (3.0 / 5.0) * (f(at + 2 * step) - f(at - 2 * step)) / (4 * step) +
           (1.0 / 10.0) * (f(at + 3 * step) - f(at - 3 * step)) / (6 * step);
}

}; // namespace diff

#endif