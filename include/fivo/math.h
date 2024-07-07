#ifndef FIVO_MATH_H
#define FIVO_MATH_H

#include <type_traits>

namespace fivo {

namespace math {

template<typename T, typename F, typename DF>
inline std::enable_if_t<std::is_arithmetic<T>::value, T>
newton_raphson(T const& x0, F&& f, DF&& df, T const& tol) {
  bool converged = false;
  T x = x0, xold = x0;
  while (!converged) {
    x = xold - f(x) / df(x);
    auto const dx = 2 * (x - xold) / (x + xold);
    if (dx < tol) converged = true;
  }
  return x;
}

} // namespace math

} // namespace fivo

#endif // FIVO_MATH_H
