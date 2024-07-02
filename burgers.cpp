#include "fivo.h"

int main(int argc, char** argv) {
  double const t0 = 0;
  double const tf = 5.0;
  double const dt = 1e-2;

  fivo::Mesh mesh(-10., 10., 500);

  // Boundary conditions
  auto left_bc = fivo::Burgers::BCWall::make();
  auto right_bc = fivo::Burgers::BCWall::make();

  // Create the system
  auto system = fivo::Burgers(mesh, left_bc, right_bc);
  using state_type = typename fivo::Burgers::state_type;

  // Initial value = gaussian
  auto const init = [&] (double const& x) { return state_type{std::exp(-0.5 * x * x)}; };
  auto X = system.create_init_state(mesh, init);

  // Output quantities
  auto const u = [&] (double const&, state_type const& s) { return s[0]; };

  // Solve and save for each numerical flux
  fivo::IOManager io("burgers_rusanov", 10, mesh);
  fivo::solve(io, system, fivo::Rusanov{}, fivo::EulerStep{}, X, t0, tf, dt, u);
  io.basename("burgers_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::HLL{}, fivo::EulerStep{}, X, t0, tf, dt, u);
}
