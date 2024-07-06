#include "fivo.h"

int main() {
  double const gamma = 1.4;
  double const t0 = 0;
  double const tf = 0.2;
  double const dt = 1e-2;

  auto mesh = fivo::Mesh(0., 1., 200);

  // Boundary conditions
  auto left_bc = fivo::system::IdealGasEuler::BCNeumann::make();
  auto right_bc = fivo::system::IdealGasEuler::BCNeumann::make();

  // Create the system
  auto system = fivo::system::IdealGasEuler(mesh, left_bc, right_bc, gamma);
  using state_type = typename fivo::system::IdealGasEuler::state_type;

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
  auto const density  = [&] (double const&, state_type const& s) { return s[0]; };
  auto const pressure = [&] (double const&, state_type const& s) { return system.pressure(s); };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const tenergy  = [&] (double const&, state_type const& s) { return s[2]; };
  auto const ienergy  = [&] (double const&, state_type const& s) {
                          auto const u = system.velocity(s);
                          return s[2] / s[0] - 0.5 * u  *u;
                        };
  auto const quantities = std::make_tuple(density, pressure, velocity, ienergy, tenergy);

  // Solve and save for each numerical flux
  auto io = fivo::IOManager("ideal_gas_euler_sod_rusanov", 1, mesh);
  auto X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::Rusanov{}, fivo::time::RK3Heun{}, X, t0, tf, dt, quantities);
  io.basename("ideal_gas_euler_sod_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLL{}, fivo::time::RK3Heun{}, X, t0, tf, dt, quantities);
  io.basename("ideal_gas_euler_sod_hllc");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLLC{}, fivo::time::RK3Heun{}, X, t0, tf, dt, quantities);
  io.basename("ideal_gas_euler_2_sod_hllc");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLLC{}, fivo::time::RK3DHeun{}, X, t0, tf, dt, quantities);
}
