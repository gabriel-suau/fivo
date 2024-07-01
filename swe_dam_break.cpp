#include "fivo.h"

int main() {
  double const grav = 9.81;
  double const t0 = 0;
  double const tf = 1;
  double const dt = 1e-3;

  fivo::Mesh mesh(0., 25., 800);

  // Boundary conditions
  auto left_bc = fivo::SWE::BCNeumann::make();
  auto right_bc = fivo::SWE::BCNeumann::make();

  // Friction model
  auto friction_model = fivo::SWE::NoFriction::make(0);

  // Create the system
  auto system = fivo::SWE(mesh, left_bc, right_bc, grav, friction_model);
  using state_type = typename fivo::SWE::state_type;

  // Initial value (h, q) as a function of space
  auto const init = [&] (double const& x) {
                      if (x < 0.5 * (mesh.xmin() + mesh.xmax())) return state_type{2., 0.};
                      else return state_type{1., 0.};
                    };

  fivo::IOManager io("swe_dam_break_rusanov", 10, mesh);
  auto X = system.create_init_state(mesh, init);
  fivo::solve(io, mesh, system, fivo::Rusanov{}, fivo::EulerStep{}, X, t0, tf, dt);
  io.basename("swe_dam_break_hll");
  X = system.create_init_state(mesh, init);
  fivo::solve(io, mesh, system, fivo::HLL{}, fivo::EulerStep{}, X, t0, tf, dt);
}
