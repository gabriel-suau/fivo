#ifndef FIVO_SOLVER_H
#define FIVO_SOLVER_H

#include <fivo/time.h>
#include <fivo/mesh.h>
#include <fivo/state.h>
#include <fivo/system.h>
#include <fivo/io.h>

namespace fivo {

template<typename System, typename Stepper, typename NumFlux, typename... Quantities>
void solve(System const& system, NumFlux const& numflux, Stepper const& stepper,
           typename System::global_state_type& X, double const& t0, double const& tf,
           double const& dt, IOManager const& io, std::tuple<Quantities...> const& quantities) {
  using global_state_type = typename System::global_state_type;

  // Init time/niter
  double t = t0;
  int niter = 0;

  auto const& left_bc_func = system.left_bc();
  auto const& right_bc_func = system.right_bc();
  auto const& mesh = system.mesh();

  // Rhs computation at each time stem
  auto const rhs =
    [&] (double const& t, global_state_type const& X) {
      global_state_type rhs_value(X.size());
      auto const left_bc = left_bc_func->compute(mesh, t, X.front(), X.back());
      auto const right_bc = right_bc_func->compute(mesh, t, X.back(), X.front());
      numflux.contribute(system, mesh, left_bc, right_bc, X, rhs_value);
      rhs_value += system.source(mesh, t, X);
      return rhs_value;
    };

  // Save the initial state
  if (io.save_frequency() > 0) io.save_state(t, niter, X, quantities);

  // Time loop
  while (t < tf) {
    // Adjust time step to ensure CFL condition is met
    auto const max_ws = system.max_wave_speed(X);
    auto step = std::min(dt, mesh.dx() / max_ws);
    // auto step = dt;
    step = ((t + step) <= tf) ? step : (tf - t);
    // Do one step
    stepper(X, t, step, rhs);
    ++niter;
    if (io.save_frequency() > 0 && niter % io.save_frequency() == 0)
      io.save_state(t, niter, X, quantities);
    t += step;
  }
}

template<typename System, typename Stepper, typename NumFlux, typename... Quantities>
void solve(IOManager const& io, System const& system, NumFlux const& numflux,
           Stepper const& stepper, typename System::global_state_type& X,
           double const& t0, double const& tf, double const& dt,
           Quantities const&... quantities) {
  return solve(io, system, numflux, stepper, X, t0, tf, dt, std::make_tuple(quantities...));
}

} // namespace fivo

#endif // FIVO_SOLVER_H
