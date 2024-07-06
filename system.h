#ifndef FIVO_SYSTEM_H
#define FIVO_SYSTEM_H

#include "state.h"
#include "mesh.h"
#include "vector.h"
#include "math.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>

namespace fivo {

namespace system {

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
    vector<state_type> out(in.size());
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
    vector<state_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (state_type const& s) { return wave_speeds(s); });
    return out;
  }
  auto max_wave_speed(global_state_type const& in) const {
    auto const ws = gwave_speeds(in);
    auto maxws = vector<value_type>(ws.size());
    auto comp = [&] (auto const& l, auto const& r) { return std::abs(l) < std::abs(r); };
    std::transform(ws.begin(), ws.end(), maxws.begin(),
                   [&] (auto const& s) {
                     return std::abs(*std::max_element(s.begin(), s.end(), comp));
                   });
    return *std::max_element(maxws.begin(), maxws.end());
  }
};

template<typename State>
struct HasAdmissible {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  virtual bool admissible(state_type const& in) const = 0;
  auto gwave_speeds(global_state_type const& in) const {
    vector<bool> out(in.size());
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
    vector<value_type> out(in.size());
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
    vector<value_type> out(in.size());
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
    vector<value_type> out(in.size());
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
  System(Mesh const& mesh,
         std::shared_ptr<typename HasBC<State>::BC> const& left_bc,
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

  using BC = System<state_type>::BC;
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

  LinearAdvection(Mesh const& mesh,
                  std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc,
                  value_type const& v)
    : System(mesh, left_bc, right_bc), m_v(v) {}

  value_type m_v;

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const&) const override { return m_v; }
  value_type velocity() const { return m_v; }
  LinearAdvection& velocity(value_type const& value) { m_v = value; return *this; }

  state_type flux(state_type const& s) const override { return m_v * s; }

  state_type wave_speeds(state_type const& /* s */) const override { return state_type{m_v}; }

  bool admissible(state_type const&) const override { return true; }
};

/* LWR MODEL FOR TRAFFIC FLOW */
struct LWRTrafficFlow : System<fivo::state<double, 1>>,
                        HasDensity<fivo::state<double, 1>>,
                        HasVelocity<fivo::state<double, 1>> {
  using state_type = fivo::state<double, 1>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using BC = System<state_type>::BC;
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

  LWRTrafficFlow(Mesh const& mesh,
                 std::shared_ptr<BC> const& left_bc,
                 std::shared_ptr<BC> const& right_bc,
                 value_type const& rmax,
                 value_type const& umax)
    : System(mesh, left_bc, right_bc), m_rmax(rmax), m_umax(umax) {}

  value_type m_rmax, m_umax;

  value_type rmax() const { return m_rmax; }
  value_type umax() const { return m_umax; }
  LWRTrafficFlow& rmax(value_type value) { m_rmax = value; return *this; }
  LWRTrafficFlow& umax(value_type value) { m_umax = value; return *this; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return m_umax * (1 - s[0] / m_rmax); }

  state_type flux(state_type const& s) const override { return s * velocity(s); }

  state_type wave_speeds(state_type const& s) const override {
    return state_type{m_umax * (1 - 2 * s[0] / m_rmax)};
  }

  bool admissible(state_type const& s) const override { return (s[0] > 0 && s[0] < m_rmax); }
};

/* BURGERS EQUATION */
struct Burgers : System<fivo::state<double, 1>>,
                 HasVelocity<fivo::state<double, 1>> {
  using state_type = fivo::state<double, 1>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using BC = System<state_type>::BC;
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

  using System::System;

  value_type velocity(state_type const& s) const override { return s[0]; }

  state_type flux(state_type const& s) const override {
    auto const& u = s[0];
    return state_type{0.5 * u * u};
  }

  state_type wave_speeds(state_type const& s) const override { return s; }

  bool admissible(state_type const&) const override { return true; }
};

/* LINEAR ACOUSTICS EQUATIONS : PRESSURE-VELOCITY FORMULATION */
struct LinearAcousticsPressure : System<fivo::state<double, 2>>,
                                 HasDensity<fivo::state<double, 2>>,
                                 HasPressure<fivo::state<double, 2>>,
                                 HasVelocity<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using BC = System<state_type>::BC;
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

  LinearAcousticsPressure(Mesh const& mesh,
                          std::shared_ptr<BC> const& left_bc,
                          std::shared_ptr<BC> const& right_bc,
                          value_type const& r0,
                          value_type const& c,
                          value_type const& u0 = 0)
    : System(mesh, left_bc, right_bc), m_r0(r0), m_c(c), m_u0(u0) {}

  value_type m_r0, m_c, m_u0;

  value_type density(state_type const& s) const override { return s[0] / (m_c * m_c); }
  value_type pressure(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1]; }

  state_type flux(state_type const& s) const override {
    return {m_u0 * s[0] + m_r0 * m_c * m_c * s[1],
            s[0] / m_r0 + m_u0 * s[1]};
  }

  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{m_u0 - m_c, m_u0 + m_c};
  }

  bool admissible(state_type const&) const override { return true; }
};

/* LINEAR ACOUSTICS EQUATIONS : DENSITY-VELOCITY FORMULATION */
struct LinearAcousticsDensity : System<fivo::state<double, 2>>,
                                HasDensity<fivo::state<double, 2>>,
                                HasPressure<fivo::state<double, 2>>,
                                HasVelocity<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  using BC = System<state_type>::BC;
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

  LinearAcousticsDensity(Mesh const& mesh,
                         std::shared_ptr<BC> const& left_bc,
                         std::shared_ptr<BC> const& right_bc,
                         value_type const& r0,
                         value_type const& c,
                         value_type const& u0 = 0)
    : System(mesh, left_bc, right_bc), m_r0(r0), m_c(c), m_u0(u0) {}

  value_type m_r0, m_c, m_u0;

  value_type density(state_type const& s) const override { return s[0]; }
  value_type pressure(state_type const& s) const override { return m_c * m_c * s[0]; }
  value_type velocity(state_type const& s) const override { return s[1]; }

  state_type flux(state_type const& s) const override {
    return {m_u0 * s[0] + m_r0 * s[1],
            m_c * m_c * s[0] + m_u0 * s[1]};
  }

  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{m_u0 - m_c, m_u0 + m_c};
  }

  bool admissible(state_type const&) const override { return true; }
};

/* SHALLOW WATER EQUATIONS */
struct SWE : System<fivo::state<double, 2>>,
             HasVelocity<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = System<state_type>::BC;
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

  SWE(Mesh const& mesh,
      std::shared_ptr<BC> const& left_bc,
      std::shared_ptr<BC> const& right_bc,
      value_type const& grav,
      std::shared_ptr<FrictionModel> const& fmodel = std::make_shared<NoFriction>(0))
    : System(mesh, left_bc, right_bc), m_grav(grav), m_fmodel(fmodel)
  {
    // Create empty topography
    topography([](value_type const&) { return 0.; },
               [](value_type const&) { return 0.; });
  }

  value_type m_grav;
  std::shared_ptr<FrictionModel> m_fmodel;
  vector<value_type> m_topo, m_gdz;

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

  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }

  state_type flux(state_type const& s) const override {
    auto const& [h, j] = s;
    return state_type{j, j * j / h + 0.5 * m_grav * h * h};
  }

  state_type wave_speeds(state_type const& s) const override {
    auto const h = s[0];
    auto const u = s[1] / h;
    auto const c = std::sqrt(m_grav * h);
    return state_type{u - c, u + c};
  }

  bool admissible(state_type const& s) const override { return s[0] > 0; }

  global_state_type source(Mesh const& mesh, value_type const& t,
                           global_state_type const& X) const override {
    global_state_type S(mesh.nx());
    for (int i = 0; i < mesh.nx(); ++i) {
      S[i] = state_type{0., - m_gdz[i] * X[i][0]} + m_fmodel->apply(m_grav, X[i]);
    }
    return S;
  }
};

/* TELEGRAPH EQUATIONS */
template<typename Sigma>
struct Telegraph : System<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  Telegraph(Mesh const& mesh,
            std::shared_ptr<BC> const& left_bc,
            std::shared_ptr<BC> const& right_bc,
            value_type const& velocity,
            Sigma const& sigma)
    : System(mesh, left_bc, right_bc), m_velocity(velocity), m_sigma(sigma)
  {}

  value_type m_velocity;
  Sigma m_sigma;

  auto velocity() const { return m_velocity; }
  auto sigma() const { return m_sigma; }
  Telegraph& velocity(value_type const& value) { m_velocity = value; return *this; }
  Telegraph& sigma(Sigma const& value) { m_sigma = value; return *this; }

  state_type flux(state_type const& s) const override {
    return state_type{m_velocity * s[0], -m_velocity * s[1]};
  }

  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{-m_velocity, m_velocity};
  }

  bool admissible(state_type const& s) const override { return (s[0] > 0 && s[1] > 0); }

  global_state_type source(Mesh const& mesh, value_type const& /* t */,
                           global_state_type const& X) const override {
    global_state_type src(mesh.nx());
    for (int i = 0; i < mesh.nx(); ++i) {
      auto const x = mesh.cell_center(i);
      auto const& s = X[i];
      src[i] = m_sigma(x, s) * state_type{s[1] - s[0], s[0] - s[1]};
    }
    return src;
  }
};

/* ISENTROPIC EULER EQUATIONS */
struct IsentropicEuler : System<fivo::state<double, 2>>,
                         HasDensity<fivo::state<double, 2>>,
                         HasVelocity<fivo::state<double, 2>>,
                         HasPressure<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = System<state_type>::BC;
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

  IsentropicEuler(Mesh const& mesh,
                  std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc,
                  value_type const& kappa,
                  value_type const& gamma)
    : System(mesh, left_bc, right_bc), m_kappa(kappa), m_gamma(gamma) {}

  value_type m_kappa, m_gamma;

  auto kappa() const { return m_kappa; }
  auto gamma() const { return m_gamma; }
  IsentropicEuler& kappa(value_type const& value) { m_kappa = value; return *this; }
  IsentropicEuler& gamma(value_type const& value) { m_gamma = value; return *this; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const override {
    return m_kappa * std::pow(s[0], m_gamma);
  }

  state_type flux(state_type const& s) const override {
    auto const p = pressure(s);
    auto const& [r, j] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p};
  }

  state_type wave_speeds(state_type const& s) const override {
    auto const p = pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_gamma * p / r);
    return state_type{u - c, u + c};
  }

  bool admissible(state_type const& s) const override { return (s[0] > 0); }
};

/* ISOTHERMAL EULER EQUATIONS */
struct IsothermalEuler : System<fivo::state<double, 2>>,
                         HasDensity<fivo::state<double, 2>>,
                         HasVelocity<fivo::state<double, 2>>,
                         HasPressure<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = System<state_type>::BC;
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

  IsothermalEuler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc, value_type const& c)
    : System(mesh, left_bc, right_bc), m_c(c) {}

  value_type m_c;

  auto c() const { return m_c; }
  IsothermalEuler& c(value_type const& value) { m_c = value; return *this; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const override { return m_c * m_c * s[0]; }

  state_type flux(state_type const& s) const override {
    auto const p = pressure(s);
    auto const& [r, j] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p};
  }

  state_type wave_speeds(state_type const& s) const override {
    auto const u = s[1] / s[0];
    return state_type{u - m_c, u + m_c};
  }

  bool admissible(state_type const& s) const override { return (s[0] > 0); }
};

/* GENERAL EULER EQUATIONS (ABSTRACT) */
struct Euler : System<fivo::state<double, 3>>,
               HasDensity<fivo::state<double, 3>>,
               HasVelocity<fivo::state<double, 3>>,
               HasPressure<fivo::state<double, 3>> {
  using state_type = fivo::state<double, 3>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = System<state_type>::BC;
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

  Euler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
        std::shared_ptr<BC> const& right_bc, value_type const& gamma)
    : System(mesh, left_bc, right_bc), m_gamma(gamma) {}

  value_type m_gamma;

  auto gamma() const { return m_gamma; }
  Euler& gamma(value_type const& value) { m_gamma = value; return *this; }

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }

  state_type flux(state_type const& s) const override {
    auto const p = pressure(s);
    auto const& [r, j, e] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p, (e + p) * u};
  }

  state_type wave_speeds(state_type const& s) const override {
    auto const p = pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_gamma * p / r);
    return state_type{u - c, u, u + c};
  }

  bool admissible(state_type const& s) const override {
    auto const r = s[0];
    auto const e = s[2];
    return (r > 0 && e > 0);
  }
};

/* GENERAL EULER EQUATIONS WITH AN IDEAL GAS EOS */
struct IdealGasEuler : Euler {
  using Euler::Euler;
  value_type pressure(state_type const& s) const override {
    auto const& [r, j, e] = s;
    return (m_gamma - 1) * (e - 0.5 * j * j / r);
  }

//   auto solve_riemann(state_type const& left, state_type const& right) const {
//     auto const& rl = left[0];
//     auto const& rr = right[0];
//     auto const pl = pressure(left);
//     auto const pr = pressure(right);
//     auto const ul = velocity(left);
//     auto const ur = velocity(right);
//     auto const cl = std::sqrt(m_gamma * pl / rl);
//     auto const cr = std::sqrt(m_gamma * pr / rr);

//     // Find p and u in the star region (pstar, ustar)
//     auto const Al = 2. / (rl * (m_gamma + 1));
//     auto const Ar = 2. / (rr * (m_gamma + 1));
//     auto const Bl = pl * (m_gamma - 1) / (m_gamma + 1);
//     auto const Br = pr * (m_gamma - 1) / (m_gamma + 1);
//     auto const du = ur - ul;

//     auto const fk =
//       [=] (auto const p, auto const pk, auto const ck, auto const Ak, auto const Bk) {
//         if (p > pk) return (p - pk) * std::sqrt(Ak / (p + Bk));
//         else return 2 * ck / (m_gamma - 1) * (std::pow(p / pk, (m_gamma - 1) / (2 * m_gamma)) - 1);
//       };
//     auto const dfk =
//       [=] (auto const p, auto const pk, auto const rk,
//           auto const ck, auto const Ak, auto const Bk) {
//         if (p > pk) return std::sqrt(Ak / (Bk + p)) * (1 - (p - pk) / (2 * (Bk + p)));
//         else return std::pow(p / pk, - (m_gamma + 1) / (2 * m_gamma)) / (rk * ck);
//       };
//     auto const fl = [&] (auto const p) { return fk(p, pl, cl, Al, Bl); };
//     auto const fr = [&] (auto const p) { return fk(p, pr, cr, Ar, Br); };
//     auto const dfl = [&] (auto const p) { return dfk(p, pl, rl, cl, Al, Bl); };
//     auto const dfr = [&] (auto const p) { return dfk(p, pr, rr, cr, Ar, Br); };
//     auto const f = [&] (auto const p) { return fl(p) + fr(p) + du; };
//     auto const df = [&] (auto const p) { return dfl(p) + dfr(p); };

//     auto const tol = 1e-6;
//     auto const ppv = 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rl + rr) * (cl + cr);
//     auto const p0 = std::max(tol, ppv);
//     auto const pstar = math::newton_raphson(p0, f, df, tol);
//     auto const ustar = 0.5 * ((ul + ur) + fr(pstar) - fl(pstar));

//     // Find rl and rr in the star region (rlstar, rrstar)
//     value_type rlstar, rrstar, sl, sr;
//     // Left shock
//     if (pstar > pl) {
//       auto const fac = (m_gamma-1)/(m_gamma+1);
//       rlstar = rl * (pstar / pl + fac) / ((pstar / pl) * fac + 1);
//       sl = ul - cl * std::sqrt((m_gamma + 1) / (2 * m_gamma) * pstar / pl + (m_gamma - 1) / (2 * m_gamma));
//     }
//     // Left rarefaction
//     else {
//       rlstar = rl * std::pow(pstar / pl, 1 / m_gamma);
//     }
//     // Right shock

//     // Right rarefaction
//   }
};

/* GENERAL EULER EQUATIONS WITH A STIFFENED GAS EOS */
struct StiffenedGasEuler : Euler {
  StiffenedGasEuler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                    std::shared_ptr<BC> const& right_bc, value_type const& gamma,
                    value_type const& p0)
    : Euler(mesh, left_bc, right_bc, gamma), m_p0(p0) {}

  value_type m_p0;
  auto p0() const { return m_p0; }
  StiffenedGasEuler& p0(value_type const& value) { m_p0 = value; return *this; }

  value_type pressure(state_type const& s) const override {
    auto const& [r, j, e] = s;
    return (m_gamma - 1) * (e - 0.5 * j * j / r) - m_gamma * m_p0;
  }
};

/* P1 RADIATIVE TRANSFER EQUATIONS */
struct RadiativeTransferP1 : System<fivo::state<double, 2>> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  RadiativeTransferP1(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
            std::shared_ptr<BC> const& right_bc,
            value_type const& c, value_type const& sigma)
    : System(mesh, left_bc, right_bc), m_c(c), m_sigma(sigma)
  {}

  value_type m_c;
  value_type m_sigma;

  auto c() const { return m_c; }
  auto sigma() const { return m_sigma; }
  RadiativeTransferP1& c(value_type const& value) { m_c = value; return *this; }
  RadiativeTransferP1& sigma(value_type const& value) { m_sigma = value; return *this; }

  state_type flux(state_type const& s) const override {
    return state_type{s[1], m_c * m_c * s[0] / 3};
  }

  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{-m_c / std::sqrt(3), m_c / std::sqrt(3)};
  }

  bool admissible(state_type const& /* s */) const override { return true; }

  global_state_type source(Mesh const& mesh, value_type const& /* t */,
                           global_state_type const& X) const override {
    global_state_type src(mesh.nx());
    std::transform(X.begin(), X.end(), src.begin(),
                   [&] (auto const& s) { return state_type{0, -m_c * m_sigma * s[1]}; });
    return src;
  }
};

/* M1 RADIATIVE TRANSFER EQUATIONS */
struct RadiativeTransferM1 : System<fivo::state<double, 3>> {
  using state_type = fivo::state<double, 3>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  RadiativeTransferM1(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
            std::shared_ptr<BC> const& right_bc,
            value_type const& c, value_type const& sigma)
    : System(mesh, left_bc, right_bc), m_c(c), m_sigma(sigma)
  {}

  value_type m_c;
  value_type m_sigma;

  auto c() const { return m_c; }
  auto sigma() const { return m_sigma; }
  RadiativeTransferM1& c(value_type const& value) { m_c = value; return *this; }
  RadiativeTransferM1& sigma(value_type const& value) { m_sigma = value; return *this; }

  state_type flux(state_type const& s) const override {
    return state_type{s[1], m_c * m_c * s[0] / 3};
  }

  state_type wave_speeds(state_type const& /* s */) const override {
    return state_type{-m_c / std::sqrt(3), m_c / std::sqrt(3)};
  }

  bool admissible(state_type const& s) const override { return (s[0] > 0 && s[1] > 0); }

  global_state_type source(Mesh const& mesh, value_type const& /* t */,
                           global_state_type const& X) const override {
    global_state_type src(mesh.nx());
    std::transform(X.begin(), X.end(), src.begin(),
                   [&] (auto const& s) { return state_type{0, -m_c * m_sigma * s[1]}; });
    return src;
  }
};

} // namespace system

} // namespace fivo

#endif // FIVO_SYSTEM_H
