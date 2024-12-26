#include "../include/fivo/fivo.h"

int main(int argc, char** argv) {
  double const t0 = 0;
  double const tf = 10.0;
  double const dt = 1e-2;

  auto mesh = fivo::Mesh(-10., 10., 500);

  // Boundary conditions
  auto left_bc = fivo::system::Burgers::BCReflective::make();
  auto right_bc = fivo::system::Burgers::BCReflective::make();

  // Create the system
  auto system = fivo::system::Burgers(mesh, left_bc, right_bc);
  using state_type = typename fivo::system::Burgers::state_type;

  // Initial value = gaussian
  auto const init = [&] (double const& x) { return state_type{std::exp(-0.5 * x * x)}; };

  // Output quantities
  auto const u = [&] (double const&, state_type const& s) { return s[0]; };
  auto const quantities = std::make_tuple(std::make_pair("velocity", u));

  // Solve and save for each numerical flux
  auto const fluxes = std::make_tuple(fivo::flux::Rusanov{},
                                      fivo::flux::HLL{},
                                      fivo::flux::Godunov{});
  fivo::traits::for_each
    (fluxes,
     [&] (auto const& flux) {
       auto io = fivo::IOManager(std::string("burgers_") + flux.name(), 10, mesh);
       auto X = fivo::discretize(mesh, init);
       fivo::solve(io, system, flux, fivo::time::RK1{}, X, t0, tf, dt, quantities);
     });
}
