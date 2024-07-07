#ifndef FIVO_MATH_H
#define FIVO_MATH_H

#include <type_traits>
#include <cmath>

namespace fivo {

namespace math {

template<typename T, typename F, typename DF>
inline std::enable_if_t<std::is_arithmetic<T>::value, T>
newton_raphson(T const& x0, F&& f, DF&& df, T const& tol) {
  bool converged = false;
  T x = x0, xold = x0;
  while (!converged) {
    x = xold - f(xold) / df(xold);
    auto const dx = 2 * std::abs(x - xold) / (x + xold);
    xold = x;
    if (dx < tol) converged = true;
  }
  return x;
}

} // namespace math

} // namespace fivo

#endif // FIVO_MATH_H
