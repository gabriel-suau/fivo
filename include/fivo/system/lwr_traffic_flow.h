#ifndef FIVO_LWR_TRAFFIC_FLOW_H
#define FIVO_LWR_TRAFFIC_FLOW_H

#include <fivo/system/system_base.h>

namespace fivo {

namespace system {

class LWRTrafficFlow : public SystemBase<fivo::state<double, 1>>,
                       public HasDensity<fivo::state<double, 1>>,
                       public HasVelocity<fivo::state<double, 1>>,
                       public HasRiemannSolver<LWRTrafficFlow, fivo::state<double, 1>> {
public:
  using base_type = SystemBase<fivo::state<double, 1>>;
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

  struct Params { value_type rmax, umax; };

  LWRTrafficFlow(Mesh const& mesh,
                 std::shared_ptr<BC> const& left_bc,
                 std::shared_ptr<BC> const& right_bc,
                 Params const& params)
    : base_type(mesh, left_bc, right_bc), m_params(params) {}

  auto const& get_params() const { return m_params; }
  void set_params(Params const& params) { m_params = params; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override {
    return m_params.umax * (1 - s[0] / m_params.rmax); }
  flux_type flux(state_type const& s) const override { return flux_type{{{{s[0] * velocity(s)}}}}; }
  state_type wave_speeds(state_type const& s) const override {
    return state_type{m_params.umax * (1 - 2 * s[0] / m_params.rmax)};
  }
  bool admissible(state_type const& s) const override {
    return (s[0] > 0 && s[0] < m_params.rmax);
  }
  state_type prim_to_cons(state_type const& s) const override { return s; }
  state_type cons_to_prim(state_type const& s) const override { return s; }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& rr = right[0];
    auto const wl = wave_speeds(left)[0];
    auto const wr = wave_speeds(right)[0];
    // Rankine-Hugoniot condition for shocks
    auto const s = m_params.umax * (1 - (rl + rr) / m_params.rmax);
    auto const fpinv = [=] (value_type const& xt) {
                         return state_type{(m_params.rmax / 2) * (1 - m_params.umax * xt)};
                       };
    // Solution of the Riemann problem
    auto sol = [=] (value_type const& xt) {
                 // Shock (concave flux)
                 if (rl < rr) { return (xt < s) ? left : right; }
                 // Rarefaction
                 if (xt < wl) return left;
                 if (xt > wr) return right;
                 return fpinv(xt);
               };
    return sol;
  }

private:
  Params m_params;
};

} // namespace system

} // namespace fivo

#endif // FIVO_LWR_TRAFFIC_FLOW_H
