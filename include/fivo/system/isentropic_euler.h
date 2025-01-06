#ifndef FIVO_ISENTROPIC_EULER_H
#define FIVO_ISENTROPIC_EULER_H

#include <fivo/system/system_base.h>

namespace fivo {

namespace system {

class IsentropicEuler : public SystemBase<fivo::state<double, 2>>,
                        public HasDensity<fivo::state<double, 2>>,
                        public HasVelocity<fivo::state<double, 2>>,
                        public HasPressure<fivo::state<double, 2>>,
                        public HasRiemannSolver<IsentropicEuler, fivo::state<double, 2>> {
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

  struct Params { value_type kappa, gamma; };

  IsentropicEuler(Mesh const& mesh,
                  std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc,
                  Params const& params)
    : base_type(mesh, left_bc, right_bc), m_params(params) {}

  auto const& get_params() const { return m_params; }
  void set_params(Params const& params) { m_params = params; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const override {
    return m_params.kappa * std::pow(s[0], m_params.gamma);
  }
  flux_type flux(state_type const& s) const override {
    auto const p = pressure(s);
    auto const& [r, j] = s;
    auto const u = j / r;
    return flux_type{{{{j}, {r * u * u + p}}}};
  }
  state_type wave_speeds(state_type const& s) const override {
    auto const p = pressure(s);
    auto const& [r, u] = cons_to_prim(s);
    auto const c = std::sqrt(m_params.gamma * p / r);
    return state_type{u - c, u + c};
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

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& rr = right[0];
    auto const ul = velocity(left);
    auto const ur = velocity(right);
    auto const& kappa = m_params.kappa;
    auto const& gamma = m_params.gamma;
    auto const cl = std::sqrt(kappa * gamma * std::pow(rl, gamma - 1));
    auto const cr = std::sqrt(kappa * gamma * std::pow(rr, gamma - 1));

    // Find hstar and ustar in the star region
    auto const du = ur - ul;
    auto const fk =
      [=] (auto const r, auto const rk) {
        if (r > rk) {
          auto const p = kappa * std::pow(r, gamma - 1);
          auto const pk = kappa * std::pow(rk, gamma - 1);
          return (r - rk) * std::sqrt((pk - p) / (r * rk * (rk - r))); 
        }
        else {
          auto const cstar = kappa * gamma * std::pow(r, gamma - 1);
          return 2 / (gamma - 1) * (cstar - cr);
        }
      };
    auto const dfk =
      [=] (auto const r, auto const rk) {
        if (r > rk) return 1.;
        else return 1.;
      };

    auto const fl = [&] (auto const r) { return fk(r, rl); };
    auto const fr = [&] (auto const r) { return fk(r, rr); };
    auto const dfl = [&] (auto const r) { return dfk(r, rl); };
    auto const dfr = [&] (auto const r) { return dfk(r, rr); };

    auto const f = [&] (auto const r) { return fl(r) + fr(r) + du; };
    auto const df = [&] (auto const r) { return dfl(r) + dfr(r); };

    auto const tol = 1e-10;
    auto const rrv = 0.5 * (rl + rr);
    auto const r0 = std::max(tol, rrv);
    auto const rstar = math::newton_raphson(r0, f, df, tol);
    auto const ustar = 0.5 * ((ul + ur) + fr(rstar) - fl(rstar));
    auto const cstar = std::sqrt(kappa * gamma * std::pow(rstar, gamma - 1));

    // Left shock speed
    auto const sl = ul - 2 / (gamma - 1) * (cstar - cl);

    // Left head and tail rarefaction speeds
    auto const shl = ul - cl;
    auto const stl = ustar - cstar;

    // Right shock speed
    auto const sr = ur - 2 / (gamma - 1) * (cstar - cr);

    // Right head and tail rarefaction speeds
    auto const shr = ur + cr;
    auto const str = ustar + cstar;

    // Solution in the star region\fan
    auto const lrstar =
      [=] (value_type const&) { return prim_to_cons(state_type{rstar, ustar}); };

    // Solution in the fan region for left rarefaction
    auto const lfan =
      [=] (value_type const& xt) {
        auto const rlfan = 1.;
        auto const ulfan = 1.;
        return prim_to_cons(state_type{rlfan, ulfan});
      };

    // Solution in the fan region for right rarefaction
    auto const rfan =
      [=] (value_type const& xt) {
        auto const rrfan = 1.;
        auto const urfan = 1.;
        return prim_to_cons(state_type{rrfan, urfan});
      };

    auto const sol =
      [=] (value_type const& xt) {
        if (xt < ustar) { // Left
          if (rstar > rl) { return (xt < sl) ? left : lrstar(xt); }
          else { return (xt < shl) ? left : (xt > stl) ? lrstar(xt) : lfan(xt); }
        }
        else {
          if (rstar > rr) { return (xt > sr) ? right : lrstar(xt); }
          else { return (xt > shr) ? right : (xt < str) ? lrstar(xt) : rfan(xt); }
        }
      };
    return sol;
  }

private:
  Params m_params;
};
                 
} // namespace system

} // namespace fivo

#endif // FIVO_ISENTROPIC_EULER_H
