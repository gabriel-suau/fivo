#include "../include/fivo/fivo.h"

int main() {
  double const gamma = 1.4;
  double const t0 = 0;
  double const tf = 0.2;
  double const dt = 1e-3;

  auto mesh = fivo::Mesh(0, 1, 400);

  // Boundary conditions
  auto left_bc = fivo::system::IdealGasEulerP<1>::BCTransmissive::make();
  auto right_bc = fivo::system::IdealGasEulerP<1>::BCTransmissive::make();

  // Create the system
  auto system = fivo::system::IdealGasEulerP<1>(mesh, left_bc, right_bc, gamma);
  using state_type = typename fivo::system::IdealGasEulerP<1>::state_type;

  // Initial value (r, q, u) as a function of space
  auto const rl = 1.0, pl = 1.0, ul = 0.0, cl = 1.;
  auto const rr = 0.125, pr = 0.1, ur = 0.0, cr = 0.125;
  auto const left = system.prim_to_cons(state_type{rl, ul, pl, cl});
  auto const right = system.prim_to_cons(state_type{rr, ur, pr, cr});
  auto const init = [&] (double const& x) {
                      return (x <= 0.5 * (mesh.xmin() + mesh.xmax())) ? left : right;
                    };

  // Output quantities
  auto const density  = [&] (double const&, state_type const& s) { return s[0]; };
  auto const pressure = [&] (double const&, state_type const& s) { return system.pressure(s); };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const tenergy  = [&] (double const&, state_type const& s) { return s[2]; };
  auto const ienergy  = [&] (double const&, state_type const& s) {
                          auto const u = system.velocity(s);
                          return s[2] / s[0] - 0.5 * u  *u;
                        };
  auto const concentration = [&] (double const&, state_type const& s) { return s[3] / s[0]; };
  auto const quantities = std::make_tuple(density, pressure, velocity, ienergy, tenergy, concentration);

  // Solve and save for each numerical flux
  auto io = fivo::IOManager("test_sod_rusanov", 1, mesh);
  auto X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::Rusanov{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("test_sod_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLL{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("test_sod_hllc");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLLC{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("test_sod_godunov");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::Godunov{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);

  // Exact solution
  auto const exact = system.solve_riemann(left, right);
  X = system.create_init_state(mesh, init);
  io.basename("test_sod_exact");
  int const nt = (tf - t0) / dt;
  for (int it = 0; it < nt; ++it) {
    auto const t = t0 + it * dt;
    for (int i = 0; i < mesh.nx(); ++i) {
      // The riemann solver is centered in 0, so we have to shift our x-position
      // when sampling the solution
      auto const x = mesh.cell_center(i) - 0.5 * (mesh.xmin() + mesh.xmax());
      X[i] = exact(x / t);
    }
    if (it % io.save_frequency() == 0) io.save_state(tf, it, X, quantities);
  }
}
