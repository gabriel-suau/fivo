#ifndef FIVO_TELEGRAPH_H
#define FIVO_TELEGRAPH_H

#include <fivo/system/system_base.h>

namespace fivo {

namespace system {

class Telegraph : public SystemBase<fivo::state<double, 2>>,
                  public HasRiemannSolver<Telegraph, fivo::state<double, 2>> {
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

  struct Params { value_type velocity, sigma; };

  Telegraph(Mesh const& mesh,
            std::shared_ptr<BC> const& left_bc,
            std::shared_ptr<BC> const& right_bc,
            Params const& params)
    : base_type(mesh, left_bc, right_bc), m_params(params) {}

  auto const& get_params() const { return m_params; }
  void set_params(Params const& params) { m_params = params; }

  flux_type flux(state_type const& s) const override {
    return flux_type{{{{m_params.velocity * s[0]}, {-m_params.velocity * s[1]}}}};
  }
  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{-m_params.velocity, m_params.velocity};
  }
  bool admissible(state_type const& s) const override { return (s[0] > 0 && s[1] > 0); }
  state_type prim_to_cons(state_type const& s) const override { return s; }
  state_type cons_to_prim(state_type const& s) const override { return s; }

  global_state_type source(Mesh const& mesh, value_type const& /* t */,
                           global_state_type const& X) const override {
    global_state_type src(mesh.nx());
    std::transform(X.begin(), X.end(), src.begin(),
                   [&] (auto const& s) {
                     return m_params.sigma * state_type{s[1] - s[0], s[0] - s[1]};
                   });
    return src;
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto exact = [&] (value_type const&) {
                   if (m_params.velocity > 0) return state_type{left[0], right[1]};
                   return state_type{left[1], right[0]};
                 };
    return exact;
  }

private:
  Params m_params;
};

} // namespace system

} // namespace fivo

#endif // FIVO_TELEGRAPH_H
