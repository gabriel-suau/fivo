#include "../include/fivo/fivo.h"

int main(int argc, char** argv) {
  double const umax = 1.;
  double const rhomax = 1.;
  double const t0 = 0;
  double const tf = 60.0;
  double const dt = 1e-2;

  auto mesh = fivo::Mesh(-10., 10., 500);

  // Boundary conditions
  auto left_bc = fivo::system::LWRTrafficFlow::BCTransmissive::make();
  auto right_bc = fivo::system::LWRTrafficFlow::BCTransmissive::make();

  // Create the system
  auto system = fivo::system::LWRTrafficFlow(mesh, left_bc, right_bc, {rhomax, umax});
  using state_type = typename fivo::system::LWRTrafficFlow::state_type;

  // Riemann problem (rarefaction)
  auto const left = state_type{0.8};
  auto const right = state_type{0.2};
  auto const init = [&] (double const& x) {
                      return (x < 0.5 * (mesh.xmin() + mesh.xmax())) ? left : right;
                    };

  // Output quantities
  auto const density = [&] (double const&, state_type const& s) { return s[0]; };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const quantities = std::make_tuple(std::make_pair("density", density),
                                          std::make_pair("velocity", velocity));

  // Solve and save for each numerical flux
  auto const fluxes = std::make_tuple(fivo::numflux::Rusanov{},
                                      fivo::numflux::HLL{},
                                      fivo::numflux::Godunov{});
  fivo::traits::for_each
    (fluxes,
     [&] (auto const& flux) {
       auto io = fivo::IOManager(std::string("traffic_") + flux.name(), 1, mesh);
       auto X = fivo::discretize(mesh, init);
       fivo::solve(io, system, flux, fivo::time::RK1{}, X, t0, tf, dt, quantities);
     });

  // Exact solution
  auto const exact = system.solve_riemann(left, right);
  auto Xe = fivo::discretize(mesh, init);
  auto io = fivo::IOManager("traffic_exact", 1, mesh);
  int const nt = (tf - t0) / dt;
  for (int it = 0; it < nt; ++it) {
    auto const t = t0 + it * dt;
    for (int i = 0; i < mesh.nx(); ++i) {
      auto const x = mesh.cell_center(i);
      Xe[i] = exact(x / t);
    }
    if (it % io.save_frequency() == 0) io.save_state(tf, it, Xe, quantities);
  }
}
