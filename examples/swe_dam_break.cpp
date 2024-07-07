#include "../include/fivo/fivo.h"

int main() {
  double const grav = 9.81;
  double const t0 = 0;
  double const tf = 1;
  double const dt = 1e-3;

  auto mesh = fivo::Mesh(0., 25., 800);

  // Boundary conditions
  auto left_bc = fivo::system::SWE::BCNeumann::make();
  auto right_bc = fivo::system::SWE::BCNeumann::make();

  // Friction model
  auto friction_model = fivo::system::SWE::NoFriction::make(0);

  // Create the system
  auto system = fivo::system::SWE(mesh, left_bc, right_bc, grav, friction_model);
  using state_type = typename fivo::system::SWE::state_type;

  // Initial value (h, q) as a function of space
  auto const init = [&] (double const& x) {
                      if (x < 0.5 * (mesh.xmin() + mesh.xmax())) return state_type{2., 0.};
                      else return state_type{1., 0.};
                    };

  // Output quantities
  auto const water_height = [&] (double const&, state_type const& s) { return s[0]; };
  auto const discharge = [&] (double const&, state_type const& s) { return s[1]; };
  auto const total_height = [&] (double const&, state_type const& s) { return s[0]; };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const topography = [&] (double const&, state_type const&) { return 0; };
  auto const quantities = std::make_tuple(topography, water_height, total_height, discharge, velocity);

  auto io = fivo::IOManager("swe_dam_break_rusanov", 10, mesh);
  auto X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::Rusanov{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("swe_dam_break_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLL{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
}
