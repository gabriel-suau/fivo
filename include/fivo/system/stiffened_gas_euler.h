#ifndef FIVO_STIFFENED_GAS_EULER_H
#define FIVO_STIFFENED_GAS_EULER_H

#include <fivo/system/euler_base.h>

namespace fivo {

namespace system {

class StiffenedGasEuler : public EulerBase {
public:
  using base_type = EulerBase;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  struct Params : base_type::Params { value_type p0; };

  StiffenedGasEuler(Mesh const& mesh,
                   std::shared_ptr<BC> const& left_bc,
                   std::shared_ptr<BC> const& right_bc,
                   Params const& params)
    : base_type(mesh, left_bc, right_bc, {}), m_params(params) {}

  auto const& get_params() const { return m_params; }
  void set_params(Params const& params) { m_params = params; }

  value_type pressure(state_type const& s) const override {
    auto const& [r, j, e] = s;
    return (m_params.gamma - 1) * (e - 0.5 * j * j / r) - m_params.gamma * m_params.p0;
  }

private:
  value_type energy(state_type const& prim) const override {
    auto const& [r, u, p] = prim;
    return (p + m_params.gamma * m_params.p0) / (m_params.gamma - 1) + 0.5 * r * u * u;
  }

  Params m_params;
};

} // namespace system

} // namespace fivo

#endif // FIVO_STIFFENED_GAS_EULER_H
