#include "fivo.h"

int main(int argc, char** argv) {
  double const u0 = 0;
  double const c = 340;
  double const r0 = 1.225;
  double const t0 = 0;
  double const tf = 0.1;

  auto mesh = fivo::Mesh(-1, 1., 200.);

  // Respect exactly the CFL condition
  auto const dt = mesh.dx() / (u0 + c);

  // Boundary conditions
  auto left_bc = fivo::system::LinearAcousticsPressure::BCNeumann::make();
  auto right_bc = fivo::system::LinearAcousticsPressure::BCNeumann::make();

  // Create the system
  auto system = fivo::system::LinearAcousticsPressure(mesh, left_bc, right_bc, c, r0, u0);
  using state_type = typename fivo::system::LinearAcousticsPressure::state_type;

  // Initial value : p = gaussian, u = 0
  auto const init = [&] (double const& x) {
                      if (x < 0.5 * (mesh.xmin() + mesh.xmax())) return state_type{1.};
                      else return state_type{0.6};
                    };
  // auto const init = [&] (double const& x) { return state_type{std::exp(-0.5 * x * x), 0.}; };
  auto X = system.create_init_state(mesh, init);

  // Output quantities
  auto const density  = [&] (double const&, state_type const& s) { return system.density(s); };
  auto const pressure = [&] (double const&, state_type const& s) { return system.pressure(s); };
  auto const velocity = [&] (double const&, state_type const& s) { return system.velocity(s); };
  auto const quantities = std::make_tuple(density, pressure, velocity);

  // Solve and save for each numerical flux
  auto io = fivo::IOManager("acoustics_rusanov", 1, mesh);
  fivo::solve(io, system, fivo::flux::Rusanov{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
  io.basename("acoustics_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, system, fivo::flux::HLL{}, fivo::time::RK1{}, X, t0, tf, dt, quantities);
}
