#ifndef FIVO_TIME_H
#define FIVO_TIME_H

#include "io.h"

namespace fivo {

struct EulerStep {
  template<typename State, typename RHS>
  void operator()(State& X, double const& t, double const& dt, RHS const& f) const {
    X += dt * f(t, X);
  }
};

struct HeunStep {
  template<typename State, typename RHS>
  void operator()(State& X, double const& t, double const& dt, RHS const& f) const {
    auto const k1 = f(t, X);
    auto const k2 = f(t + dt, X + dt * k1);
    X += dt * 0.5 * (k1 + k2);
  }
};

struct MidPointStep {
  template<typename State, typename RHS>
  void operator()(State& X, double const& t, double const& dt, RHS const& f) const {
    auto const k1 = f(t, X);
    auto const k2 = f(t + 0.5 * dt, X + 0.5 * dt * k1);
    X += dt * k2;
  }
};

template<typename StepImpl, typename State, typename RHS>
void time_solve(IOManager const& io, StepImpl const& step_impl, State& X,
                double const& t0, double const& tf,
                double const& dt, RHS const& f) {
  double t = t0;
  int niter = 0;
  while (t < tf) {
    auto const step = ((t + dt) <= tf) ? dt : (tf - t);
    step_impl(X, t, step, f);
    ++niter;
    if (niter % io.save_frequency() == 0) io.save_state(t, X);
    t += step;
  }
}

} // namespace fivo

#endif // FIVO_TIME_H
