#ifndef FIVO_LINEAR_ACOUSTICS_H
#define FIVO_LINEAR_ACOUSTICS_H

#include <fivo/system/system_base.h>

#include <fivo/riemann.h>

namespace fivo {

namespace system {

class LinearAcousticsPressure : public SystemBase<fivo::state<double, 2>>,
                                public HasDensity<fivo::state<double, 2>>,
                                public HasPressure<fivo::state<double, 2>>,
                                public HasVelocity<fivo::state<double, 2>>,
                                public HasRiemannSolver<LinearAcousticsPressure,
                                                        fivo::state<double, 2>> {
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

  struct Params { value_type c, r0, u0; };

  LinearAcousticsPressure(Mesh const& mesh,
                          std::shared_ptr<BC> const& left_bc,
                          std::shared_ptr<BC> const& right_bc,
                          Params const& params)
    : base_type(mesh, left_bc, right_bc), m_c(params.c), m_r0(params.r0), m_u0(params.u0) {}

  auto get_params() const { return Params{m_c, m_r0, m_u0}; }
  void set_params(Params const& params) {
    m_c = params.c;
    m_r0 = params.r0;
    m_u0 = params.u0;
  }

  value_type density(state_type const& s) const override {
    return s[0] / (m_c * m_c);
  }
  value_type pressure(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1]; }
  flux_type flux(state_type const& s) const override {
    return flux_type{{{{m_u0 * s[0] + m_r0 * m_c * m_c * s[1]}, {s[0] / m_r0 + m_u0 * s[1]}}}};
  }
  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{m_u0 - m_c, m_u0 + m_c};
  }
  bool admissible(state_type const&) const override { return true; }
  state_type prim_to_cons(state_type const& s) const override { return s; }
  state_type cons_to_prim(state_type const& s) const override { return s; }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& pl = left[0];
    auto const& ul = left[1];
    auto const& pr = right[0];
    auto const& ur = right[1];
    auto const pstar = 0.5 * (pl + pr - (ur - ul) * (m_r0 * m_c));
    auto const ustar = 0.5 * (ul + ur - (pr - pl) / (m_r0 * m_c));
    auto const wr = m_u0 - m_c;
    auto const wl = m_u0 + m_c;
    // Solution of the Riemann problem
    auto sol = [=] (value_type const& xt) {
                 if (xt < wl) return left;
                 if (xt > wr) return right;
                 return state_type{pstar, ustar};
               };
    return sol;
  }

  auto wave_directions(state_type const& s) const {
    auto const z0 = m_r0 * m_c;
    return array<state_type, 2>{state_type{-z0, 1}, state_type{z0, 1}};
  }

  friend struct ::fivo::riemann::Exact;

private:
  value_type m_c, m_r0, m_u0;
};

} // namespace system

template<>
inline auto
riemann::Exact::solve(system::LinearAcousticsPressure const& sys,
                      typename system::LinearAcousticsPressure::state_type const& left,
                      typename system::LinearAcousticsPressure::state_type const& right,
                      array<double, 1> const& normal) {
  using state_type = typename system::LinearAcousticsPressure::state_type;
  using value_type = typename system::LinearAcousticsPressure::value_type;
  auto const z0 = sys.m_r0 * sys.m_c;
  auto const speeds = sys.wave_speeds(left);
  auto const directions = sys.wave_directions(right);
  auto const delta = right - left;
  auto const a1 = (-delta[0] + z0 * delta[1]) / (2 * z0);
  auto const a2 = (delta[0] + z0 * delta[1]) / (2 * z0);
  auto const waves = array<state_type, 2> {a1 * directions[0], a2 * directions[1]};
  auto eval = [&] (value_type const& xt) {
                if (xt < speeds[0]) { return left; }
                if (xt > speeds[1]) { return right; }
                return state_type{};
              };
  return std::make_tuple(waves, speeds, eval);
}

namespace system {

/* LINEAR ACOUSTICS EQUATIONS : DENSITY-VELOCITY FORMULATION */
class LinearAcousticsDensity : public SystemBase<fivo::state<double, 2>>,
                               public HasDensity<fivo::state<double, 2>>,
                               public HasPressure<fivo::state<double, 2>>,
                               public HasVelocity<fivo::state<double, 2>>,
                               public HasRiemannSolver<LinearAcousticsDensity,
                                                       fivo::state<double, 2>> {
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

  struct Params { value_type c, r0, u0; };

  LinearAcousticsDensity(Mesh const& mesh,
                         std::shared_ptr<BC> const& left_bc,
                         std::shared_ptr<BC> const& right_bc,
                         Params const& params)
    : base_type(mesh, left_bc, right_bc), m_params(params) {}

  auto const& get_params() const { return m_params; }
  void set_params(Params const& params) { m_params = params; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type pressure(state_type const& s) const override {
    return m_params.c * m_params.c * s[0];
  }
  value_type velocity(state_type const& s) const override { return s[1]; }
  flux_type flux(state_type const& s) const override {
    auto const& [c, r0, u0] = m_params;
    return flux_type{{{{u0 * s[0] + r0 * s[1]},
            {c * c / r0 * s[0] + u0 * s[1]}}}};
  }
  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{m_params.u0 - m_params.c, m_params.u0 + m_params.c};
  }
  bool admissible(state_type const&) const override { return true; }
  state_type prim_to_cons(state_type const& s) const override { return s; }
  state_type cons_to_prim(state_type const& s) const override { return s; }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& ul = left[1];
    auto const& rr = right[0];
    auto const& ur = right[1];
    auto const& [c, r0, u0] = m_params;
    auto const rstar = 0.5 * (rl + rr - (ur - ul) * r0 / c);
    auto const ustar = 0.5 * (rl + ur - (rr - rl) * c / r0);
    auto const wr = u0 - c;
    auto const wl = u0 + c;
    // Solution of the Riemann problem
    auto sol = [=] (value_type const& xt) {
                 if (xt < wl) return left;
                 if (xt > wr) return right;
                 return state_type{rstar, ustar};
               };
    return sol;
  }

private:
  Params m_params;
};

} // namespace system

} // namespace fivo

#endif // FIVO_LINEAR_ACOUSTICS_H
