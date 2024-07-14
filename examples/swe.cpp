#include "../include/fivo/fivo.h"

int main() {
  double const grav = 9.81;
  double const t0 = 0;
  double const tf = 40;
  double const dt = 2e-3;

  auto mesh = fivo::Mesh(0., 25., 400);

  // Boundary conditions
  // auto left_bc_func = [] (double const& t) { return state_type{2. + std::sin(t / (2. * M_PI)), 1.}; };
  // auto left_bc = fivo::SWE::BCImposedWaterHeight<decltype(left_bc_func)>::make(left_bc_func);
  auto left_bc = fivo::system::SWE::BCReflective::make();
  auto right_bc = fivo::system::SWE::BCReflective::make();

  // Friction model
  auto friction_model = fivo::system::SWE::ManningFriction::make(0.04);
  // auto friction_model = fivo::SWE::DarcyWeisbachFriction::make(0.093);

  // Create the system
  auto system = fivo::system::SWE(mesh, left_bc, right_bc, grav, friction_model);
  using state_type = typename fivo::system::SWE::state_type;

  // Topography as a function of space
  auto const bump = [] (double const& x) {
                      if (8. < x && x < 12.) return 0.2 - 0.05 * (x - 10) *  (x - 10);
                      return 0.;
                    };
  auto const dbump_dx = [] (double const& x) {
                          if (8. < x && x < 12.) return - 0.05 * 2. * (x - 10);
                          return 0.;
                        };
  system.topography(bump, dbump_dx);

  // Initial value (h, q) as a function of space
  auto const init = [&] (double const& x) { return state_type{std::max(0.3 + 0.02*std::sin(x) - bump(x), 0.), 0.}; };
  auto X = system.create_init_state(mesh, init);

  // Output quantities
  auto const topography = [&] (double const& x, state_type const&) { return bump(x); };
  auto const water_height = [&] (double const&, state_type const& s) { return s[0]; };
  auto const total_height = [&] (double const& x, state_type const& s) { return s[0] + bump(x); };
  auto const discharge = [&] (double const&, state_type const& s) { return s[1]; };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const quantities = std::make_tuple(std::make_pair("topography", topography),
                                          std::make_pair("water_height", water_height),
                                          std::make_pair("total_height", total_height),
                                          std::make_pair("discharge", discharge),
                                          std::make_pair("velocity", velocity));

  // Solve and save for each numerical flux
  auto io = fivo::IOManager("swe_rusanov", 20, mesh);
  io.basename("swe_rusanov");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::Rusanov{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("swe_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLL{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
}
