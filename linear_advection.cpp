#include "fivo.h"

int main(int argc, char** argv) {
  double const v = 1.;
  double const t0 = 0;
  double const tf = 10.0;
  double const dt = 1e-3;

  fivo::Mesh mesh(-5., 5., 1000);

  // Boundary conditions
  auto left_bc = fivo::LinearAdvection::BCPeriodic::make();
  auto right_bc = fivo::LinearAdvection::BCPeriodic::make();

  // Create the system
  auto system = fivo::LinearAdvection(mesh, left_bc, right_bc, v);
  using state_type = typename fivo::LinearAdvection::state_type;

  // Initial value = gaussian
  auto const init = [&] (double const& x) { return state_type{std::exp(-0.5 * x * x)}; };
  auto X = system.create_init_state(mesh, init);

  // Output quantities
  auto const u = [&] (double const&, state_type const& s) { return s[0]; };

  // Solve and save for each numerical flux
  fivo::IOManager io("advection_rusanov", 10, mesh);
  fivo::solve(io, system, fivo::Rusanov{}, fivo::HeunStep{}, X, t0, tf, dt, u);
  io.basename("advection_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::HLL{}, fivo::HeunStep{}, X, t0, tf, dt, u);
}
