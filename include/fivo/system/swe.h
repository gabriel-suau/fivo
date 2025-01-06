#ifndef FIVO_SWE_H
#define FIVO_SWE_H

#include <fivo/system/system_base.h>

namespace fivo {

namespace system {

class SWE : public SystemBase<fivo::state<double, 2>>,
            public HasVelocity<fivo::state<double, 2>>,
            public HasRiemannSolver<SWE, fivo::state<double, 2>> {
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

  template<typename Func>
  struct BCImposedWaterHeight : BC {
    Func func;
    static inline auto make(Func const& func) { return std::make_shared<BCImposedWaterHeight>(func); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      return func(t);
    }
  };
  template<typename Func>
  struct BCImposedDischarge : BC {
    Func func;
    static inline auto make(Func const& func) { return std::make_shared<BCImposedDischarge>(func); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      return func(t);
    }
  };

  /* FRICTION MODELS */
  struct FrictionModel {
    explicit FrictionModel(value_type c) : m_c(c) {}
    value_type m_c;
    virtual state_type apply(value_type const& /* g */,
                             state_type const& /* s */) const = 0;
  };
  struct NoFriction : FrictionModel {
    using FrictionModel::FrictionModel;
    state_type apply(value_type const& g,
                     state_type const& s) const final {
      return {0., 0.};
    }
    static inline auto make(value_type c) { return std::make_shared<NoFriction>(c); }
  };
  struct LaminarFriction : FrictionModel {
    using FrictionModel::FrictionModel;
    state_type apply(value_type const& g,
                     state_type const& s) const final {
      auto const& h = s[0];
      auto const u = s[1] / s[0];
      return {0., m_c * u / (g * h * h)};
    }
    static inline auto make(value_type c) { return std::make_shared<LaminarFriction>(c); }
  };
  struct ManningFriction : FrictionModel {
    using FrictionModel::FrictionModel;
    state_type apply(value_type const& /* g */,
                     state_type const& s) const final {
      auto const& h = s[0];
      auto const u = s[1] / s[0];
      return {0., m_c * m_c * std::abs(u) / std::pow(h, 4. / 3.)};
    }
    static inline auto make(value_type c) { return std::make_shared<ManningFriction>(c); }
  };
  struct DarcyWeisbachFriction : FrictionModel {
    using FrictionModel::FrictionModel;
    state_type apply(value_type const& g,
                     state_type const& s) const final {
      auto const& h = s[0];
      auto const u = s[1] / s[0];
      return {0., m_c * std::abs(u) / (8 * g * h)};
    }
    static inline auto make(value_type c) { return std::make_shared<DarcyWeisbachFriction>(c); }
  };

  struct Params { value_type grav; };

  SWE(Mesh const& mesh,
      std::shared_ptr<BC> const& left_bc,
      std::shared_ptr<BC> const& right_bc,
      Params const& params,
      std::shared_ptr<FrictionModel> const& fmodel = std::make_shared<NoFriction>(0))
    : base_type(mesh, left_bc, right_bc), m_params(params), m_fmodel(fmodel)
  {
    // Create empty topography
    topography([](value_type const&) { return 0.; },
               [](value_type const&) { return 0.; });
  }

  auto const& get_params() const { return m_params; }
  void set_params(Params const& params) { m_params = params; }

  auto friction_model() const { return m_fmodel; }
  SWE& friction_model(std::shared_ptr<FrictionModel> const& value) { m_fmodel = value; return *this; }

  auto topography() const { return m_topo; }
  template<typename TopoFunc, typename TopoDxFunc>
  SWE& topography(TopoFunc&& func, TopoDxFunc&& dfunc) {
    auto const nx = m_mesh.nx();
    auto const dx = m_mesh.dx();
    // Create topo and gdz
    m_topo.resize(nx);
    m_gdz.resize(nx);
    for (int i = 0; i < nx; ++i) {
      auto const x = m_mesh.xmin() + (i + 0.5) * dx;
      m_topo[i] = func(x);
      m_gdz[i] = m_params.grav * dfunc(x);
    }
    return *this;
  }
  template<typename TopoFunc>
  SWE& topography(TopoFunc&& func) {
    auto const nx = m_mesh.nx();
    auto const dx = m_mesh.dx();
    // Create topo
    m_topo.resize(nx);
    for (int i = 0; i < nx; ++i) {
      auto const x = m_mesh.xmin() + (i + 0.5) * dx;
      m_topo[i] = func(x);
    }
    // Create gdz
    m_gdz.resize(nx);
    m_gdz[0] = m_params.grav * (-3 * m_topo[0] + 4 * m_topo[1] - m_topo[2]) / (2 * dx);
    for (int i = 1; i < nx - 1; ++i) {
      m_gdz[i] = m_params.grav * (m_topo[i+1] - m_topo[i-1]) / (2 * dx);
    }
    m_gdz[nx - 1] = m_params.grav * (3 * m_topo[nx-1] - 4 * m_topo[nx-2] + m_topo[nx-3]) / (2 * dx);
    return *this;
  }

  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }
  flux_type flux(state_type const& s) const override {
    auto const& [h, j] = s;
    return flux_type{{{{j}, {j * j / h + 0.5 * m_params.grav * h * h}}}};
  }
  state_type wave_speeds(state_type const& s) const override {
    auto const h = s[0];
    auto const u = s[1] / h;
    auto const c = std::sqrt(m_params.grav * h);
    return state_type{u - c, u + c};
  }
  bool admissible(state_type const& s) const override { return s[0] > 0; }
  state_type prim_to_cons(state_type const& s) const override {
    auto const& [h, u] = s;
    return state_type{h, h * u};
  }
  state_type cons_to_prim(state_type const& s) const override {
    auto const& [h, j] = s;
    return state_type{h, j / h};
  }

  global_state_type source(Mesh const& mesh, value_type const& t,
                           global_state_type const& X) const override {
    global_state_type S(mesh.nx());
    for (int i = 0; i < mesh.nx(); ++i) {
      S[i] = state_type{0., - m_gdz[i] * X[i][0]} + m_fmodel->apply(m_params.grav, X[i]);
    }
    return S;
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& g = m_params.grav;
    auto const& hl = left[0];
    auto const& hr = right[0];
    auto const ul = velocity(left);
    auto const ur = velocity(right);
    auto const cl = std::sqrt(g * hl);
    auto const cr = std::sqrt(g * hr);

    // Find hstar and ustar in the star region
    auto const du = ur - ul;
    auto const fv =
      [=] (auto const h, auto const hk) { return std::sqrt(0.5 * g * (1 / h + 1 / hk)); };

    auto const fk =
      [=] (auto const h, auto const hk) {
        if (h > hk) return (h - hk) * fv(h, hk);
        else return 2 * (std::sqrt(g * h) - std::sqrt(g * hk));
      };
    auto const dfk =
      [=] (auto const h, auto const hk) {
        if (h > hk) return fv(h, hk) + g * (h - hk) * (hk - h) / (2 * h * h * fv(h, hk));
        else return std::sqrt(g / h);
      };

    auto const fl = [&] (auto const h) { return fk(h, hl); };
    auto const fr = [&] (auto const h) { return fk(h, hr); };
    auto const dfl = [&] (auto const h) { return dfk(h, hl); };
    auto const dfr = [&] (auto const h) { return dfk(h, hr); };

    auto const f = [&] (auto const h) { return fl(h) + fr(h) + du; };
    auto const df = [&] (auto const h) { return dfl(h) + dfr(h); };

    auto const tol = 1e-10;
    auto const hhv = 0.5 * (hl + hr);
    auto const h0 = std::max(tol, hhv);
    auto const hstar = math::newton_raphson(h0, f, df, tol);
    auto const ustar = 0.5 * ((ul + ur) + fr(hstar) - fl(hstar));
    auto const cstar = std::sqrt(g * hstar);

    // Left shock speed
    auto const sl = ul - std::sqrt(0.5 * g * hstar * hl * (hstar + hl)) / hl;

    // Left head and tail rarefaction speeds
    auto const shl = ul - cl;
    auto const stl = ustar - cstar;

    // Right shock speed
    auto const sr = ur + std::sqrt(0.5 * g * hstar * hr * (hstar + hr)) / hr;

    // Right head and tail rarefaction speeds
    auto const shr = ur + cr;
    auto const str = ustar + cstar;

    // Solution in the star region\fan
    auto const lrstar =
      [=] (value_type const&) { return prim_to_cons(state_type{hstar, ustar}); };

    // Solution in the fan region for left rarefaction
    auto const lfan =
      [=] (value_type const& xt) {
        auto const hlfan = std::pow(2 * cl + ul - xt, 2) / (9 * g);
        auto const ulfan = ul + 2. / 3. * (xt - ul + cl);
        return prim_to_cons(state_type{hlfan, ulfan});
      };

    // Solution in the fan region for right rarefaction
    auto const rfan =
      [=] (value_type const& xt) {
        auto const hrfan = std::pow(2 * cr - ur + xt, 2) / (9 * g);
        auto const urfan = ur + 2. / 3. * (xt - ur - cr);
        return prim_to_cons(state_type{hrfan, urfan});
      };

    auto const sol =
      [=] (value_type const& xt) {
        if (xt < ustar) { // Left
          if (hstar > hl) { return (xt < sl) ? left : lrstar(xt); }
          else { return (xt < shl) ? left : (xt > stl) ? lrstar(xt) : lfan(xt); }
        }
        else {
          if (hstar > hr) { return (xt > sr) ? right : lrstar(xt); }
          else { return (xt > shr) ? right : (xt < str) ? lrstar(xt) : rfan(xt); }
        }
      };
    return sol;
  }

private:
  Params m_params;
  std::shared_ptr<FrictionModel> m_fmodel;
  vector<value_type> m_topo, m_gdz;
};

} // namespace system

} // namespace fivo

#endif // FIVO_SWE_H
