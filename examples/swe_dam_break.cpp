#include "../include/fivo/fivo.h"

int main() {
  double const grav = 9.81;
  double const t0 = 0;
  double const tf = 1;
  double const dt = 1e-3;

  auto mesh = fivo::Mesh(0., 25., 200);

  // Boundary conditions
  auto left_bc = fivo::system::SWE::BCTransmissive::make();
  auto right_bc = fivo::system::SWE::BCTransmissive::make();

  // Friction model
  auto friction_model = fivo::system::SWE::NoFriction::make(0);

  // Create the system
  auto system = fivo::system::SWE(mesh, left_bc, right_bc, grav, friction_model);
  using state_type = typename fivo::system::SWE::state_type;

  // Initial value (h, q) as a function of space
  auto const left = state_type{2, 0};
  auto const right = state_type{1, 0};
  auto const init = [&] (double const& x) {
                      if (x < 0.5 * (mesh.xmin() + mesh.xmax())) return left;
                      else return right;
                    };

  // Output quantities
  auto const water_height = [&] (double const&, state_type const& s) { return s[0]; };
  auto const discharge = [&] (double const&, state_type const& s) { return s[1]; };
  auto const total_height = [&] (double const&, state_type const& s) { return s[0]; };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const topography = [&] (double const&, state_type const&) { return 0; };
  auto const quantities = std::make_tuple(std::make_pair("topography", topography),
                                          std::make_pair("water_height", water_height),
                                          std::make_pair("total_height", total_height),
                                          std::make_pair("discharge", discharge),
                                          std::make_pair("velocity", velocity));

  // Solve and save for each numerical flux
  auto const fluxes = std::make_tuple(fivo::flux::Rusanov{},
                                      fivo::flux::HLL{},
                                      fivo::flux::Godunov{});
  fivo::traits::for_each
    (fluxes,
     [&] (auto const& flux) {
       auto io = fivo::IOManager(std::string("dam_break_") + flux.name(), 1, mesh);
       auto X = fivo::discretize(mesh, init);
       fivo::solve(io, system, flux, fivo::time::RK1{}, X, t0, tf, dt, quantities);
     });

  // Exact solution
  auto const exact = system.solve_riemann(left, right);
  auto X = fivo::discretize(mesh, init);
  auto io = fivo::IOManager("dam_break_exact", 1, mesh);
  int const nt = (tf - t0) / dt;
  for (int it = 0; it < nt; ++it) {
    auto const t = t0 + it * dt;
    for (int i = 0; i < mesh.nx(); ++i) {
      // The riemann solver is centered in 0, so we have to shift our x-position
      // when sampling the solution
      auto const x = mesh.cell_center(i) - 0.5 * (mesh.xmin() + mesh.xmax());
      X[i] = exact(x / t);
    }
    if (it % io.save_frequency() == 0) io.save_state(t, it+1, X, quantities);
  }
}
