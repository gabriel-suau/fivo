#include "../include/fivo/fivo.h"

int main(int argc, char** argv) {
  double const v = 1.;
  double const t0 = 0;
  double const tf = 10.0;

  auto mesh = fivo::Mesh(-5., 5., 100);

  // Exactly meet the CFL condition
  double const dt = mesh.dx() / v;

  // Boundary conditions
  auto left_bc = fivo::system::LinearAdvection::BCPeriodic::make();
  auto right_bc = fivo::system::LinearAdvection::BCPeriodic::make();

  // Create the system
  auto system = fivo::system::LinearAdvection(mesh, left_bc, right_bc, v);
  using state_type = typename fivo::system::LinearAdvection::state_type;

  // Initial value = gaussian
  auto const init = [&] (double const& x) { return state_type{std::exp(-0.5 * x * x)}; };
  auto X = system.create_init_state(mesh, init);

  // Output quantities
  auto const u = [&] (double const&, state_type const& s) { return s[0]; };

  // Solve and save for each numerical flux
  auto io = fivo::IOManager("advection_rusanov", 1, mesh);
  fivo::solve(io, system, fivo::flux::Rusanov{}, fivo::time::RK1{}, X, t0, tf, dt, u);
  io.basename("advection_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLL{}, fivo::time::RK1{}, X, t0, tf, dt, u);
  io.basename("advection_godunov");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::Godunov{}, fivo::time::RK1{}, X, t0, tf, dt, u);
}
