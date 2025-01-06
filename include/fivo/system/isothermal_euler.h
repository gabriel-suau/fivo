#ifndef FIVO_ISOTHERMAL_EULER_H
#define FIVO_ISOTHERMAL_EULER_H

#include <fivo/system/system_base.h>

namespace fivo {

namespace system {

class IsothermalEuler : public SystemBase<fivo::state<double, 2>>,
                        public HasDensity<fivo::state<double, 2>>,
                        public HasVelocity<fivo::state<double, 2>>,
                        public HasPressure<fivo::state<double, 2>> {
public:
  using base_type = SystemBase<fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
  using flux_type = typename base_type::flux_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = typename base_type::BC;
  struct BCReflective : BC {
    static inline auto make() { return std::make_shared<BCReflective>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      return state_type{in[0], -in[1]};
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

  struct Params { value_type c; };

  IsothermalEuler(Mesh const& mesh,
                  std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc,
                  Params const& params)
    : base_type(mesh, left_bc, right_bc), m_params(params) {}

  auto const& get_params() const { return m_params; }
  void set_params(Params const& params) { m_params = params; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const override {
    return m_params.c * m_params.c * s[0];
  }
  flux_type flux(state_type const& s) const override {
    auto const p = pressure(s);
    auto const& [r, j] = s;
    auto const u = j / r;
    return flux_type{{{{j}, {r * u * u + p}}}};
  }
  state_type wave_speeds(state_type const& s) const override {
    auto const u = s[1] / s[0];
    return state_type{u - m_params.c, u + m_params.c};
  }
  bool admissible(state_type const& s) const override { return (s[0] > 0); }
  state_type prim_to_cons(state_type const& s) const override {
    auto const& [r, u] = s;
    return state_type{r, r * u};
  }
  state_type cons_to_prim(state_type const& s) const override {
    auto const& [r, j] = s;
    return state_type{r, j / r};
  }

private:
  Params m_params;
};

} // namespace system

} // namespace fivo

#endif // FIVO_ISOTHERMAL_EULER_H
