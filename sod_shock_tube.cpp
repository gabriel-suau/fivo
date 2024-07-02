#include "fivo.h"

int main() {
  double const gamma = 1.4;
  double const t0 = 0;
  double const tf = 0.2;
  double const dt = 1e-4;

  fivo::Mesh mesh(0., 1., 400);

  // Boundary conditions
  auto left_bc = fivo::IdealGasEuler::BCNeumann::make();
  auto right_bc = fivo::IdealGasEuler::BCNeumann::make();

  // Create the system
  auto system = fivo::IdealGasEuler(mesh, left_bc, right_bc, gamma);
  using state_type = typename fivo::IdealGasEuler::state_type;

  // Initial value (r, q, u) as a function of space
  auto const rl = 1.0, pl = 1.0, ul = 0.0;
  auto const rr = 0.125, pr = 0.1, ur = 0.0;
  auto const jl = rl * ul;
  auto const jr = rr * ur;
  auto const el = pl / (gamma - 1) + 0.5 * jl * jl * rl;
  auto const er = pr / (gamma - 1) + 0.5 * jr * jr * rr;
  auto const init = [&] (double const& x) {
                      if (x <= 0.5 * (mesh.xmin() + mesh.xmax())) return state_type{rl, jl, el};
                      else return state_type{rr, jr, er};
                    };

  // Output quantities
  auto const rho = [&] (double const&, state_type const& s) { return s[0]; };
  auto const pressure = [&] (double const&, state_type const& s) { return system.pressure(s); };
  auto const velocity = [&] (double const&, state_type const& s) { return s[1] / s[0]; };
  auto const quantities = std::make_tuple(rho, pressure, velocity);

  fivo::IOManager io("euler_sod_rusanov", 1, mesh);
  auto X = system.create_init_state(mesh, init);
  fivo::solve(io, mesh, system, fivo::Rusanov{}, fivo::EulerStep{}, X, t0, tf, dt, quantities);
  io.basename("euler_sod_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, mesh, system, fivo::HLL{}, fivo::EulerStep{}, X, t0, tf, dt, quantities);
}
