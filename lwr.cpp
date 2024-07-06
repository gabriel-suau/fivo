#include "fivo.h"

int main(int argc, char** argv) {
  double const umax = 1.;
  double const rhomax = 1.;
  double const t0 = 0;
  double const tf = 60.0;
  double const dt = 1e-2;

  auto mesh = fivo::Mesh(-10., 10., 500);

  // Boundary conditions
  auto left_bc = fivo::system::LWRTrafficFlow::BCWall::make();
  auto right_bc = fivo::system::LWRTrafficFlow::BCNeumann::make();

  // Create the system
  auto system = fivo::system::LWRTrafficFlow(mesh, left_bc, right_bc, rhomax, umax);
  using state_type = typename fivo::system::LWRTrafficFlow::state_type;

  // Initial value = gaussian
  // auto const init = [&] (double const& x) { return state_type{0.9 + 0.1 * std::exp(-0.5 * x * x)}; };
  auto const init = [&] (double const& x) { return state_type{0.9 - 0.1 * std::exp(-0.5 * x * x)}; };
  // auto const init = [&] (double const& x) { return state_type{1.}; };
  auto X = system.create_init_state(mesh, init);

  // Output quantities
  auto const density = [&] (double const&, state_type const& s) { return s[0]; };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const quantities = std::make_tuple(density, velocity);

  // Solve and save for each numerical flux
  auto io = fivo::IOManager("traffic_rusanov", 5, mesh);
  fivo::solve(io, system, fivo::flux::Rusanov{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("traffic_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLL{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("traffic_godunov");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::Godunov{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
}
