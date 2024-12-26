#ifndef FIVO_MATH_H
#define FIVO_MATH_H

#include <fivo/state.h>

#include <type_traits>
#include <cmath>
#include <algorithm>

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

template<typename State>
inline auto error_l1(global_state<State> const& X,
                     global_state<State> const& Y,
                     double const& dx) {
  State err = {0.};
  for (std::size_t i = 0; i < X.size(); ++i)
    for (std::size_t v = 0; v < X[i].size(); ++v)
      err[v] += std::abs(X[i][v] - Y[i][v]);
  return err * dx;
}

template<typename State>
inline auto error_l2(global_state<State> const& X,
                     global_state<State> const& Y,
                     double const& dx) {
  State err = {0.};
  for (std::size_t i = 0; i < X.size(); ++i) {
    auto const diff = X[i] - Y[i];
    for (std::size_t v = 0; v < X[i].size(); ++v) {
      err[v] += diff[v] * diff[v];
    }
  }
  std::for_each(err.begin(), err.end(), [](auto& e) { e = std::sqrt(e); });
  return err * std::sqrt(dx);
}

template<typename State>
inline auto error_linf(global_state<State> const& X,
                       global_state<State> const& Y,
                       double const& /* dx */) {
  State err = {0.};
  for (std::size_t i = 0; i < X.size(); ++i)
    for (std::size_t v = 0; v < X[i].size(); ++v)
      err[v] = std::max(err[v], std::abs(X[i][v] - Y[i][v]));
  return err;
}

} // namespace math

} // namespace fivo

#endif // FIVO_MATH_H
