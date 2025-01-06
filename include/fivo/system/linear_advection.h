#ifndef FIVO_LINEAR_ADVECTION_H
#define FIVO_LINEAR_ADVECTION_H

#include <fivo/system/system_base.h>

#include <fivo/riemann.h>

namespace fivo {

namespace system {

class LinearAdvection : public SystemBase<fivo::state<double, 1>>,
                        public HasDensity<fivo::state<double, 1>>,
                        public HasVelocity<fivo::state<double, 1>>,
                        public HasRiemannSolver<LinearAdvection, fivo::state<double, 1>> {
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

  struct Params { value_type velocity; };

  LinearAdvection(Mesh const& mesh,
                  std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc,
                  Params const& params)
    : base_type(mesh, left_bc, right_bc), m_velocity(params.velocity) {}

  auto get_params() const { return Params{m_velocity}; }
  void set_params(Params const& params) { m_velocity = params.velocity; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const&) const  override{ return m_velocity; }
  flux_type flux(state_type const& s) const override { return flux_type{{{{m_velocity * s[0]}}}}; }
  bool admissible(state_type const&) const override { return true; }
  state_type prim_to_cons(state_type const& s) const override { return s; }
  state_type cons_to_prim(state_type const& s) const override { return s; }

  state_type wave_speeds(state_type const& /* s */) const override { return {m_velocity}; }
  auto wave_directions(state_type const& /* s */) const {
    return array<state_type, 1>{state_type{1}};
  }
  auto wave_structure(state_type const& s) const {
    return std::make_pair(wave_directions(s), wave_speeds(s));
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    return [&] (value_type const& xt) { return (m_velocity > xt) ? left : right; };
  }

  // Exact riemann solver must be able to access private members directly
  friend struct ::fivo::riemann::Exact;

private:
  value_type m_velocity;
};

} // namespace system

/**
 * Specializes fivo::riemann::Exact::solve for the LinearAdvection system
 */
template<>
inline auto
riemann::Exact::solve(system::LinearAdvection const& sys,
                      typename system::LinearAdvection::state_type const& left,
                      typename system::LinearAdvection::state_type const& right,
                      array<double, 1> const& /* /normal */) {
  using state_type = typename system::LinearAdvection::state_type;
  using value_type = typename system::LinearAdvection::value_type;
  auto const speeds = array<value_type, 1>{sys.m_velocity};
  auto const waves = array<state_type, 1>{right - left};
  auto const apdq = (sys.m_velocity > 0) ? (speeds[0] * waves[0]) : state_type{0};
  auto const amdq = (sys.m_velocity < 0) ? (speeds[0] * waves[0]) : state_type{0};
  auto eval = [&] (value_type const& xt) { return (sys.m_velocity > xt) ? left : right; };
  return std::make_tuple(waves, speeds, apdq, amdq, eval);
}

} // namespace fivo

#endif // FIVO_LINEAR_ADVECTION_H
