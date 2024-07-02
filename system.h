#ifndef FIVO_SYSTEM_H
#define FIVO_SYSTEM_H

#include "state.h"
#include "mesh.h"
#include <cmath>
#include <algorithm>
#include <memory>

namespace fivo {

struct HasMesh {
  HasMesh() = default;
  explicit HasMesh(Mesh const& mesh) : m_mesh(mesh) {}
  Mesh const& mesh() const { return m_mesh; }
  Mesh m_mesh;
};

template<typename State>
struct HasFlux {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual state_type flux(state_type const& in) const = 0;
  auto gflux(global_state_type const& in) const {
    std::vector<state_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (state_type const& s) { return flux(s); });
    return out;
  }
};

template<typename State>
struct HasWaveSpeed {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual state_type wave_speeds(state_type const& in) const = 0;
  auto gwave_speeds(global_state_type const& in) const {
    std::vector<state_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (state_type const& s) { return wave_speeds(s); });
    return out;
  }
};

template<typename State>
struct HasAdmissible {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual bool admissible(state_type const& in) const = 0;
  auto gwave_speeds(global_state_type const& in) const {
    std::vector<bool> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (state_type const& s) { return admissible(s); });
    return out;
  }
};

template<typename State>
struct HasVelocity {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual value_type velocity(state_type const& in) const = 0;
  auto gvelocity(global_state_type const& in) const {
    std::vector<value_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (state_type const& s) { return velocity(s); });
    return out;
  }
};

template<typename State>
struct HasDensity {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual value_type density(state_type const& in) const = 0;
  auto gdensity(global_state_type const& in) const {
    std::vector<value_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (state_type const& s) { return density(s); });
    return out;
  }
};

template<typename State>
struct HasPressure {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual value_type pressure(state_type const& in) const = 0;
  auto gpressure(global_state_type const& in) const {
    std::vector<value_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (state_type const& s) { return pressure(s); });
    return out;
  }
};

template<typename State>
struct HasInitState {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  template<typename InitFunc>
  auto create_init_state(Mesh const& mesh, InitFunc&& func) const {
    global_state_type X0(mesh.nx());
    for (int i = 0; i < mesh.nx(); ++i) { X0[i] = func(mesh.cell_center(i)); }
    return X0;
  }
};

template<typename State>
struct HasSource {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual global_state_type source(Mesh const& mesh, value_type const& t,
                                   global_state_type const& X) const {
    return global_state_type(mesh.nx(), {0});
  }
};

template<typename State>
struct HasBC {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  struct BC {
    using state_type = State;
    using global_state_type = fivo::global_state<state_type>;
    using value_type = typename state_type::value_type;
    virtual state_type compute(Mesh const& mesh, value_type const& t,
                               state_type const& in, state_type const& in_opbound,
                               int const normal) const = 0;
  };

  HasBC(std::shared_ptr<BC> left_bc, std::shared_ptr<BC> right_bc)
    : m_left_bc(left_bc), m_right_bc(right_bc)
  {}

  auto left_bc() const { return m_left_bc; }
  auto right_bc() const { return m_right_bc; }
  std::shared_ptr<BC> m_left_bc, m_right_bc;
};

template<typename State>
struct System : HasMesh, HasFlux<State>, HasWaveSpeed<State>, HasAdmissible<State>, HasInitState<State>, HasSource<State>, HasBC<State> {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;
  System(Mesh const& mesh, std::shared_ptr<typename HasBC<State>::BC> const& left_bc,
         std::shared_ptr<typename HasBC<State>::BC> const& right_bc)
    : HasMesh(mesh), HasBC<State>(left_bc, right_bc) {}
};

/* LINEAR ADVECTION EQUATION */
struct LinearAdvection : System<fivo::state<double, 1>>,
                         HasDensity<fivo::state<double, 1>>,
                         HasVelocity<fivo::state<double, 1>> {
  using state_type = fivo::state<double, 1>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return -in;
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };

  LinearAdvection(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc, value_type const& v)
    : System(mesh, left_bc, right_bc), m_v(v) {}

  value_type m_v;

  value_type density(state_type const& s) const final { return s[0]; }
  value_type velocity(state_type const&) const final { return m_v; }
  value_type velocity() const { return m_v; }
  LinearAdvection& velocity(value_type const& value) { m_v = value; return *this; }

  state_type flux(state_type const& s) const final { return m_v * s; }

  state_type wave_speeds(state_type const& /* s */) const final { return state_type{m_v}; }

  bool admissible(state_type const&) const final { return true; }
};

/* LINEAR ACOUSTICS EQUATIONS : PRESSURE-VELOCITY FORMULATION */
struct LinearAcousticsPressure : System<fivo::state<double, 2>>,
                                 HasDensity<fivo::state<double, 2>>,
                                 HasPressure<fivo::state<double, 2>>,
                                 HasVelocity<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return state_type{in[0], -in[1]};
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };

  LinearAcousticsPressure(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                          std::shared_ptr<BC> const& right_bc,
                          value_type const& r0, value_type const& c, value_type const& u0 = 0)
    : System(mesh, left_bc, right_bc), m_r0(r0), m_c(c), m_u0(u0) {}

  value_type m_r0, m_c, m_u0;

  value_type density(state_type const& s) const final { return s[0] / (m_c * m_c); }
  value_type pressure(state_type const& s) const final { return s[0]; }
  value_type velocity(state_type const& s) const final { return s[1]; }

  state_type flux(state_type const& s) const final {
    return {m_u0 * s[0] + m_r0 * m_c * m_c * s[1],
            s[0] / m_r0 + m_u0 * s[1]};
  }

  state_type wave_speeds(state_type const& /* s */) const final {
    return state_type{m_u0 - m_c, m_u0 + m_c};
  }

  bool admissible(state_type const&) const final { return true; }
};

/* LINEAR ACOUSTICS EQUATIONS : DENSITY-VELOCITY FORMULATION */
struct LinearAcousticsDensity : System<fivo::state<double, 2>>,
                                HasDensity<fivo::state<double, 2>>,
                                HasPressure<fivo::state<double, 2>>,
                                HasVelocity<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return state_type{in[0], -in[1]};
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };

  LinearAcousticsDensity(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                          std::shared_ptr<BC> const& right_bc,
                          value_type const& r0, value_type const& c, value_type const& u0 = 0)
    : System(mesh, left_bc, right_bc), m_r0(r0), m_c(c), m_u0(u0) {}

  value_type m_r0, m_c, m_u0;

  value_type density(state_type const& s) const final { return s[0]; }
  value_type pressure(state_type const& s) const final { return m_c * m_c * s[0]; }
  value_type velocity(state_type const& s) const final { return s[1]; }

  state_type flux(state_type const& s) const final {
    return {m_u0 * s[0] + m_r0 * s[1],
            m_c * m_c * s[0] + m_u0 * s[1]};
  }

  state_type wave_speeds(state_type const& /* s */) const final {
    return state_type{m_u0 - m_c, m_u0 + m_c};
  }

  bool admissible(state_type const&) const final { return true; }
};

/* BURGERS EQUATION */
struct Burgers : System<fivo::state<double, 1>>,
                 HasVelocity<fivo::state<double, 1>> {
  using state_type = fivo::state<double, 1>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using System::System;

  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return -in;
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };

  value_type velocity(state_type const& s) const final { return s[0]; }

  state_type flux(state_type const& s) const final {
    auto const& u = s[0];
    return state_type{0.5 * u * u};
  }

  state_type wave_speeds(state_type const& s) const final { return s; }

  bool admissible(state_type const&) const final { return true; }
};

/* SHALLOW WATER EQUATIONS */
struct SWE : System<fivo::state<double, 2>>,
             HasVelocity<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return state_type{in[0], -in[1]};
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };
  template<typename Func>
  struct BCImposedWaterHeight : BC {
    Func func;
    static inline auto make(Func const& func) { return std::make_shared<BCImposedWaterHeight>(func); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return func(t);
    }
  };
  template<typename Func>
  struct BCImposedDischarge : BC {
    Func func;
    static inline auto make(Func const& func) { return std::make_shared<BCImposedDischarge>(func); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
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

  SWE(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
      std::shared_ptr<BC> const& right_bc,
      value_type const& grav, std::shared_ptr<FrictionModel> const& fmodel)
    : System(mesh, left_bc, right_bc), m_grav(grav), m_fmodel(fmodel)
  {
    // Create empty topography
    topography([](value_type const&) { return 0.; },
               [](value_type const&) { return 0.; });
  }

  value_type m_grav;
  std::shared_ptr<FrictionModel> m_fmodel;
  std::vector<value_type> m_topo, m_gdz;

  auto grav() const { return m_grav; }
  SWE& grav(value_type const& value) { m_grav = value; return *this; }

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
      m_gdz[i] = m_grav * dfunc(x);
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
    m_gdz[0] = m_grav * (-3 * m_topo[0] + 4 * m_topo[1] - m_topo[2]) / (2 * dx);
    for (int i = 1; i < nx - 1; ++i) {
      m_gdz[i] = m_grav * (m_topo[i+1] - m_topo[i-1]) / (2 * dx);
    }
    m_gdz[nx - 1] = m_grav * (3 * m_topo[nx-1] - 4 * m_topo[nx-2] + m_topo[nx-3]) / (2 * dx);
    return *this;
  }

  value_type velocity(state_type const& s) const final { return s[1] / s[0]; }

  state_type flux(state_type const& s) const override {
    auto const& [h, j] = s;
    return state_type{j, j * j / h + 0.5 * m_grav * h * h};
  }

  state_type wave_speeds(state_type const& s) const final {
    auto const h = s[0];
    auto const u = s[1] / h;
    auto const c = std::sqrt(m_grav * h);
    return state_type{u - c, u + c};
  }

  bool admissible(state_type const& s) const final { return s[0] > 0; }

  global_state_type source(Mesh const& mesh, value_type const& t, global_state_type const& X) const final {
    global_state_type S(mesh.nx());
    for (int i = 0; i < mesh.nx(); ++i) {
      S[i] = state_type{0., - m_gdz[i] * X[i][0]} + m_fmodel->apply(m_grav, X[i]);
    }
    return S;
  }
};

struct IsentropicEuler : System<fivo::state<double, 2>>,
                         HasDensity<fivo::state<double, 2>>,
                         HasVelocity<fivo::state<double, 2>>,
                         HasPressure<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return state_type{in[0], -in[1]};
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };

  explicit IsentropicEuler(Mesh const& mesh,
                           std::shared_ptr<BC> const& left_bc,
                           std::shared_ptr<BC> const& right_bc,
                           value_type const& kappa, value_type const& gamma)
    : System(mesh, left_bc, right_bc), m_kappa(kappa), m_gamma(gamma) {}

  value_type m_kappa, m_gamma;

  auto kappa() const { return m_kappa; }
  auto gamma() const { return m_gamma; }
  IsentropicEuler& kappa(value_type const& value) { m_kappa = value; return *this; }
  IsentropicEuler& gamma(value_type const& value) { m_gamma = value; return *this; }

  value_type density(state_type const& s) const final { return s[0]; }
  value_type velocity(state_type const& s) const final { return s[1] / s[0]; }

  value_type pressure(state_type const& s) const final {
    return m_kappa * std::pow(s[0], m_gamma);
  }

  state_type flux(state_type const& s) const final {
    auto const p = pressure(s);
    auto const& [r, j] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p};
  }

  state_type wave_speeds(state_type const& s) const final {
    auto const p = pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_gamma * p / r);
    return state_type{u - c, u + c};
  }

  bool admissible(state_type const& s) const final {
    auto const r = s[0];
    auto const p = pressure(s);
    return (r > 0 && p > 0);
  }
};

struct IdealGasEuler : System<fivo::state<double, 3>>,
                       HasDensity<fivo::state<double, 3>>,
                       HasVelocity<fivo::state<double, 3>>,
                       HasPressure<fivo::state<double, 3>> {
  using state_type = fivo::state<double, 3>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return state_type{in[0], -in[1], in[2]};
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };

  explicit IdealGasEuler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                         std::shared_ptr<BC> const& right_bc, value_type const& gamma)
    : System(mesh, left_bc, right_bc), m_gamma(gamma) {}

  value_type m_gamma;

  auto gamma() const { return m_gamma; }
  IdealGasEuler& gamma(value_type const& value) { m_gamma = value; return *this; }

  value_type density(state_type const& s) const final { return s[0]; }
  value_type velocity(state_type const& s) const final { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const final {
    auto const& [r, j, e] = s;
    return (m_gamma - 1) * (e - 0.5 * j * j / r);
  }

  state_type flux(state_type const& s) const final {
    auto const p = pressure(s);
    auto const& [r, j, e] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p, (e + p) * u};
  }

  state_type wave_speeds(state_type const& s) const final {
    auto const p = pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_gamma * p / r);
    return state_type{u - c, u, u + c};
  }

  bool admissible(state_type const& s) const final {
    auto const r = s[0];
    auto const p = pressure(s);
    return (r > 0 && p > 0);
  }
};

struct StiffenedGasEuler : System<fivo::state<double, 3>>,
                           HasDensity<fivo::state<double, 3>>,
                           HasVelocity<fivo::state<double, 3>>,
                           HasPressure<fivo::state<double, 3>> {
  using state_type = fivo::state<double, 3>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = HasBC<state_type>::BC;
  struct BCWall : BC {
    static inline auto make() { return std::make_shared<BCWall>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return state_type{in[0], -in[1], in[2]};
    }
  };
  struct BCNeumann : BC {
    static inline auto make() { return std::make_shared<BCNeumann>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in;
    }
  };
  struct BCPeriodic : BC {
    static inline auto make() { return std::make_shared<BCPeriodic>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound,
                       int const normal) const final {
      return in_opbound;
    }
  };

  explicit StiffenedGasEuler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                             std::shared_ptr<BC> const& right_bc,
                             value_type const& gamma, value_type const& p0)
    : System(mesh, left_bc, right_bc), m_gamma(gamma), m_p0(p0) {}

  value_type m_gamma;
  value_type m_p0;

  auto gamma() const { return m_gamma; }
  auto p0() const { return m_p0; }
  StiffenedGasEuler& gamma(value_type const& value) { m_gamma = value; return *this; }
  StiffenedGasEuler& p0(value_type const& value) { m_p0 = value; return *this; }

  value_type density(state_type const& s) const final { return s[0]; }
  value_type velocity(state_type const& s) const final { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const final {
    auto const& [r, j, e] = s;
    return (m_gamma - 1) * (e - 0.5 * j * j / r) - m_gamma * m_p0;
  }

  state_type flux(state_type const& s) const final {
    auto const p = pressure(s);
    auto const& [r, j, e] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p, (e + p) * u};
  }

  state_type wave_speeds(state_type const& s) const final {
    auto const p = pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_gamma * p / r);
    return state_type{u - c, u, u + c};
  }

  bool admissible(state_type const& s) const final {
    auto const r = s[0];
    auto const p = pressure(s);
    return (r > 0 && p > 0);
  }
};

} // namespace fivo

#endif // FIVO_SYSTEM_H
