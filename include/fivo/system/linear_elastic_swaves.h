#ifndef FIVO_LINEAR_ELASTIC_SWAVES_H
#define FIVO_LINEAR_ELASTIC_SWAVES_H

#include <fivo/system/system_base.h>

namespace fivo {

namespace system {

/* ELASTIC COMPRESSION WAVES IN A SOLID */
class LinearElasticSwaves : public SystemBase<fivo::state<double, 2>> {
public:
  using base_type = SystemBase<fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
  using flux_type = typename base_type::flux_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  using BC = typename base_type::BC;
  struct BCReflective : BC {
    static inline auto make() { return std::make_shared<BCReflective>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      return -in;
    }
  };
  struct BCTransmissive : BC {
    static inline auto make() { return std::make_shared<BCTransmissive>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      return in_opbound;
    }
  };

  struct Params { value_type rho, mu; };

  LinearElasticSwaves(Mesh const& mesh,
                      std::shared_ptr<BC> const& left_bc,
                      std::shared_ptr<BC> const& right_bc,
                      Params const& params)
    : base_type(mesh, left_bc, right_bc), m_rho(params.rho), m_mu(params.mu) {}

  auto get_params() const { return Params{m_rho, m_mu}; }
  void set_params(Params const& params) {
    m_rho = params.rho;
    m_mu = params.mu;
  }

  flux_type flux(state_type const& s) const override {
    return flux_type{{{{- 0.5 * s[1]}, {- 2 * m_mu / m_rho * s[0]}}}};
  }
  state_type wave_speeds(state_type const& /* s */) const override {
    auto const c = std::sqrt(m_mu / m_rho);
    return state_type{-c, c};
  }
  bool admissible(state_type const&) const override { return true; }
  state_type prim_to_cons(state_type const& s) const override { return s; }
  state_type cons_to_prim(state_type const& s) const override { return s; }

private:
  value_type m_rho, m_mu;
};

} // namespace system

} // namespace fivo

#endif // FIVO_LINEAR_ELASTIC_SWAVES_H
