#ifndef FIVO_BURGERS_H
#define FIVO_BURGERS_H

#include <fivo/system/system_base.h>

#include <fivo/riemann.h>

namespace fivo {

namespace system {

class Burgers : public SystemBase<fivo::state<double, 1>>,
                public HasVelocity<fivo::state<double, 1>>,
                public HasRiemannSolver<Burgers, fivo::state<double, 1>> {
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

  using SystemBase::SystemBase;

  value_type velocity(state_type const& s) const override { return s[0]; }
  flux_type flux(state_type const& s) const override {
    auto const& u = s[0];
    return flux_type{{{{0.5 * u * u}}}};
  }
  state_type wave_speeds(state_type const& s) const override { return s; }
  bool admissible(state_type const&) const override { return true; }
  state_type prim_to_cons(state_type const& s) const override { return s; }
  state_type cons_to_prim(state_type const& s) const override { return s; }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& ul = left[0];
    auto const& ur = right[0];
    auto const wl = wave_speeds(left)[0];
    auto const wr = wave_speeds(right)[0];
    // Rankine-Hugoniot condition for shocks
    auto const s = 0.5 * (ul + ur);
    // (f')^{-1}(xt)
    auto const fpinv = [=] (value_type const& xt) { return state_type{xt}; };
    // Solution of the Riemann problem
    auto sol = [=] (value_type const& xt) {
                 // Shock (convex flux)
                 if (ul > ur) { return (xt < s) ? left : right; }
                 // Rarefaction
                 if (xt < wl) return left;
                 if (xt > wr) return right;
                 return fpinv(xt);
               };
    return sol;
  }
};

} // namespace system

template<>
inline auto
riemann::Exact::solve(system::Burgers const& sys,
                      typename system::Burgers::state_type const& left,
                      typename system::Burgers::state_type const& right,
                      array<double, 1> const& normal) {
  using state_type = typename system::Burgers::state_type;
  using value_type = typename system::Burgers::value_type;
  auto const& ul = left[0];
  auto const& ur = right[0];
  auto const wl = sys.wave_speeds(left)[0];
  auto const wr = sys.wave_speeds(right)[0];
  auto const speeds = array<value_type, 1>{0.5 * (ul + ul)};
  auto const waves = array<state_type, 1>{right - left};
  auto const apdq = (speeds[0] > 0) ? (speeds[0] * waves[0]) : state_type{0};
  auto const amdq = (speeds[0] < 0) ? (speeds[0] * waves[0]) : state_type{0};
  auto eval = [&] (value_type const& xt) {
                if (ul > ur) { return (xt < speeds[0]) ? left : right; }
                if (xt < wl) { return left; }
                if (xt > wr) { return right; }
                return state_type{xt};
              };
  return std::make_tuple(waves, speeds, apdq, amdq, eval);
}

} // namespace fivo

#endif // FIVO_BURGERS_H
