#ifndef FIVO_SOLVER_H
#define FIVO_SOLVER_H

#include "time.h"
#include "mesh.h"
#include "state.h"
#include "system.h"
#include "io.h"

namespace fivo {

template<typename System, typename StepImpl, typename NumFlux, typename... Quantities>
void solve(IOManager const& io, System const& system, NumFlux const& numflux,
           StepImpl const& step_impl, typename System::global_state_type& X,
           double const& t0, double const& tf, double const& dt,
           std::tuple<Quantities...> const& quantities) {
  using global_state_type = typename System::global_state_type;

  // Init time/niter
  double t = t0;
  int niter = 0;

  auto const& left_bc_func = system.left_bc();
  auto const& right_bc_func = system.right_bc();
  auto const& mesh = system.mesh();

  // Rhs computation at each time stem
  auto const rhs = [&] (double const& t, global_state_type const& X) {
                     global_state_type rhs_value(X.size());
                     // Add flux contrib
                     auto const left_bc = left_bc_func->compute(mesh, t, X.front(), X.back(), 1);
                     auto const right_bc = right_bc_func->compute(mesh, t, X.back(), X.front(), -1);
                     auto const nf = numflux.gcompute(system, left_bc, right_bc, X);
                     rhs_value[0] = rhs_value[0] + nf[0];
                     for (int i = 1; i < X.size(); ++i) {
                       rhs_value[i - 1] -= nf[i];
                       rhs_value[i] += nf[i];
                     }
                     rhs_value.back() -= nf.back();
                     // Divide by mesh step
                     rhs_value /= mesh.dx();
                     // Source contrib
                     rhs_value += system.source(mesh, t, X);
                     return rhs_value;
                   };

  while (t < tf) {
    auto const step = ((t + dt) <= tf) ? dt : (tf - t);
    step_impl(X, t, step, rhs);
    ++niter;
    if (niter % io.save_frequency() == 0) io.save_state(t, niter, X, quantities);
    t += step;
  }
}

template<typename System, typename StepImpl, typename NumFlux, typename... Quantities>
void solve(IOManager const& io, System const& system, NumFlux const& numflux,
           StepImpl const& step_impl, typename System::global_state_type& X,
           double const& t0, double const& tf, double const& dt,
           Quantities const&... quantities) {
  return solve(io, system, numflux, step_impl, X, t0, tf, dt, std::make_tuple(quantities...));
}

template<typename System, typename StepImpl, typename NumFlux, typename... Quantities>
void solve_splitting(IOManager const& io, System const& system, NumFlux const& numflux,
                     StepImpl const& step_impl, typename System::global_state_type& X,
                     double const& t0, double const& tf, double const& dt,
                     std::tuple<Quantities...> const& quantities) {
  using global_state_type = typename System::global_state_type;

  // Init time/niter
  double t = t0;
  int niter = 0;

  auto const& left_bc_func = system.left_bc();
  auto const& right_bc_func = system.right_bc();
  auto const& mesh = system.mesh();

  // Flux rhs
  auto const rhs_flux =
    [&] (double const& t, global_state_type const& X) {
      global_state_type rhs_value(X.size());
      // Add flux contrib
      auto const left_bc = left_bc_func->compute(mesh, t, X.front(), X.back(), 1);
      auto const right_bc = right_bc_func->compute(mesh, t, X.back(), X.front(), -1);
      auto const nf = numflux.gcompute(system, left_bc, right_bc, X);
      rhs_value[0] = rhs_value[0] + nf[0];
      for (int i = 1; i < X.size(); ++i) {
        rhs_value[i - 1] -= nf[i];
        rhs_value[i] += nf[i];
      }
      rhs_value.back() -= nf.back();
      // Divide by mesh step
      rhs_value /= mesh.dx();
      return rhs_value;
    };

  // Source rhs
  auto const rhs_src = [&] (double const& t, global_state_type const& X) {
                         return system.source(mesh, t, X);
                       };

  // Save initial solution
  io.save_state(t, niter, X, quantities);

  // Time loop
  while (t < tf) {
    auto const step = ((t + dt) <= tf) ? dt : (tf - t);
    step_impl(X, t, step, rhs_flux);
    step_impl(X, t, step, rhs_src);
    ++niter;
    if (niter % io.save_frequency() == 0) io.save_state(t, niter, X, quantities);
    t += step;
  }
}

} // namespace fivo

#endif // FIVO_SOLVER_H
