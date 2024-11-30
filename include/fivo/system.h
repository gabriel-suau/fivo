#ifndef FIVO_SYSTEM_H
#define FIVO_SYSTEM_H

#include "state.h"
#include "mesh.h"
#include "vector.h"
#include "math.h"
#include "traits.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>

namespace fivo {

namespace system {

template<typename System, typename State>
struct HasVelocity {
  using state_type = State;
  using value_type = typename state_type::value_type;

  value_type velocity(state_type const& in) const {
    return static_cast<System const*>(this)->velocity(in);
  }
};

template<typename System, typename State>
struct HasDensity {
  using state_type = State;
  using value_type = typename state_type::value_type;

  value_type density(state_type const& in) const {
    return static_cast<System const*>(this)->density(in);
  }
};

template<typename System, typename State>
struct HasPressure {
  using state_type = State;
  using value_type = typename state_type::value_type;

  value_type pressure(state_type const& in) const {
    return static_cast<System const*>(this)->pressure(in);
  }
};

template<typename System, typename State>
struct HasRiemannSolver {
  using state_type = State;
  using value_type = typename state_type::value_type;

  auto solve_riemann(state_type const& left, state_type const& right) const {
    return static_cast<System const*>(this)->solve_riemann(left, right);
  }
};

} // namespace system

namespace traits {

template<typename System>
struct has_velocity {
  static constexpr bool value =
    traits::is_derived_v<System, system::HasVelocity<System, typename System::state_type>>;
};
template<typename System>
constexpr bool has_velocity_v = has_velocity<System>::value;

template<typename System>
struct has_density {
  static constexpr bool value =
    traits::is_derived_v<System, system::HasDensity<System, typename System::state_type>>;
};
template<typename System>
constexpr bool has_density_v = has_density<System>::value;

template<typename System>
struct has_pressure {
  static constexpr bool value =
    traits::is_derived_v<System, system::HasPressure<System, typename System::state_type>>;
};
template<typename System>
constexpr bool has_pressure_v = has_pressure<System>::value;

template<typename System>
struct has_riemann_solver {
  static constexpr bool value =
    traits::is_derived_v<System, system::HasRiemannSolver<System, typename System::state_type>>;
};
template<typename System>
constexpr bool has_riemann_solver_v = has_riemann_solver<System>::value;

} // namespace traits

namespace system {

template<typename Derived, typename State>
struct System {
  using state_type = State;
  using global_state_type = fivo::global_state<State>;
  using value_type = typename state_type::value_type;

  /* BOUNDARY CONDITIONS BASE CLASS */
  struct BC {
    using state_type = State;
    using global_state_type = fivo::global_state<state_type>;
    using value_type = typename state_type::value_type;
    virtual state_type compute(Mesh const& mesh, value_type const& t,
                               state_type const& in, state_type const& in_opbound) const = 0;
  };

  System(Mesh const& mesh,
         std::shared_ptr<BC> const& left_bc,
         std::shared_ptr<BC> const& right_bc)
    : m_mesh(mesh), m_left_bc(left_bc), m_right_bc(right_bc) {}

  /* LOCAL STATE FUNCTIONS */
  /**
   * \brief Computes the pysical flux for a given local state.
   *
   * \param[in] in The input local state
   * \return The physical flux
   */
  state_type flux(state_type const& in) const {
    return static_cast<Derived const*>(this)->flux(in);
  }

  /**
   * \brief Computes the wave speeds (i.e. the eigenvalues of the flux jacobian)
   * for a given local state.
   *
   * \param[in] in The input local state
   * \return All wave speeds
   */
  state_type wave_speeds(state_type const& in) const {
    return static_cast<Derived const*>(this)->wave_speeds(in);
  }

  /**
   * \brief Determines if a local state is admissible
   *
   * \param[in] in The input local state
   * \return True or False depending on the admissibility conditions
   */
  bool admissible(state_type const& in) const {
    return static_cast<Derived const*>(this)->admissible(in);
  }

  /**
   * \brief Transforms a local state in conservative variables to primitives variables.
   *
   * For some systems, this may be the identity function.
   *
   * \param[in] in The input local state in conservative variables
   * \return The same physical local state in primitive variables formulation
   */
  state_type cons_to_prim(state_type const& in) const {
    return static_cast<Derived const*>(this)->cons_to_prim(in);
  }

  /**
   * \brief Transforms a local state in primitive variables to conservative variables.
   *
   * For some systems, this may be the identity function.
   *
   * \param[in] in The input local state in primitive variables
   * \return The same physical local state in conservative variables formulation
   */
  state_type prim_to_cons(state_type const& in) const {
    return static_cast<Derived const*>(this)->prim_to_cons(in);
  }

  /* GLOBAL STATE FUNCTION */
  template<typename LocalFunc>
  auto global(global_state_type const& in, LocalFunc&& func) const {
    using ret_type = decltype(func(std::declval<state_type>()));
    vector<ret_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(), func);
    return out;
  }

  auto max_wave_speed(global_state_type const& in) const {
    auto comp = [&] (auto const& l, auto const& r) { return std::abs(l) < std::abs(r); };
    value_type max = 0;
    for (auto const& s : in) {
      auto const ws = wave_speeds(s);
      for (auto const& w : ws) {
        auto const aw = std::abs(w);
        if (aw > max) max = aw;
      }
    }
    return max;
  }

  template<typename InitFunc>
  auto create_init_state(Mesh const& mesh, InitFunc&& func) const {
    global_state_type X0(mesh.nx());
    for (int i = 0; i < mesh.nx(); ++i) { X0[i] = func(mesh.cell_center(i)); }
    return X0;
  }

  virtual global_state_type source(Mesh const& mesh, value_type const& t,
                                   global_state_type const& X) const {
    return global_state_type(mesh.nx(), {0});
  }

  auto const& mesh() const { return m_mesh; }
  auto const& left_bc() const { return m_left_bc; }
  auto const& right_bc() const { return m_right_bc; }

protected:
  Mesh m_mesh;
  std::shared_ptr<BC> m_left_bc, m_right_bc;
};

/* LINEAR ADVECTION EQUATION */
struct LinearAdvection : System<LinearAdvection, fivo::state<double, 1>>,
                         HasDensity<LinearAdvection, fivo::state<double, 1>>,
                         HasVelocity<LinearAdvection, fivo::state<double, 1>>,
                         HasRiemannSolver<LinearAdvection, fivo::state<double, 1>> {
  using base_type = System<LinearAdvection, fivo::state<double, 1>>;
  using state_type = typename base_type::state_type;
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

  LinearAdvection(Mesh const& mesh,
                  std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc,
                  value_type const& v)
    : base_type(mesh, left_bc, right_bc), m_v(v) {}

  value_type m_v;

  value_type density(state_type const& s) const { return s[0]; }
  value_type velocity(state_type const&) const { return m_v; }
  value_type velocity() const { return m_v; }
  LinearAdvection& velocity(value_type const& value) { m_v = value; return *this; }

  state_type flux(state_type const& s) const { return m_v * s; }
  state_type wave_speeds(state_type const& /* s */) const { return state_type{m_v}; }
  bool admissible(state_type const&) const { return true; }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    return [&] (value_type const&) { return  (m_v > 0) ? left : right; };
  }
};

/* LWR MODEL FOR TRAFFIC FLOW */
struct LWRTrafficFlow : System<LWRTrafficFlow, fivo::state<double, 1>>,
                        HasDensity<LWRTrafficFlow, fivo::state<double, 1>>,
                        HasVelocity<LWRTrafficFlow, fivo::state<double, 1>>,
                        HasRiemannSolver<LWRTrafficFlow, fivo::state<double, 1>> {
  using base_type = System<LWRTrafficFlow, fivo::state<double, 1>>;
  using state_type = typename base_type::state_type;
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

  LWRTrafficFlow(Mesh const& mesh,
                 std::shared_ptr<BC> const& left_bc,
                 std::shared_ptr<BC> const& right_bc,
                 value_type const& rmax,
                 value_type const& umax)
    : base_type(mesh, left_bc, right_bc), m_rmax(rmax), m_umax(umax) {}

  value_type m_rmax, m_umax;

  value_type rmax() const { return m_rmax; }
  value_type umax() const { return m_umax; }

  value_type density(state_type const& s) const { return s[0]; }
  value_type velocity(state_type const& s) const { return m_umax * (1 - s[0] / m_rmax); }
  state_type flux(state_type const& s) const { return s * velocity(s); }
  state_type wave_speeds(state_type const& s) const {
    return state_type{m_umax * (1 - 2 * s[0] / m_rmax)};
  }
  bool admissible(state_type const& s) const { return (s[0] > 0 && s[0] < m_rmax); }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& rr = right[0];
    auto const wl = wave_speeds(left)[0];
    auto const wr = wave_speeds(right)[0];
    // Rankine-Hugoniot condition for shocks
    auto const s = m_umax * (1 - (rl + rr) / m_rmax);
    auto const fpinv =
      [=] (value_type const& xt) { return state_type{(m_rmax / 2) * (1 - m_umax * xt)}; };
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
};

/* BURGERS EQUATION */
struct Burgers : System<Burgers, fivo::state<double, 1>>,
                 HasVelocity<Burgers, fivo::state<double, 1>>,
                 HasRiemannSolver<Burgers, fivo::state<double, 1>> {
  using base_type = System<Burgers, fivo::state<double, 1>>;
  using state_type = typename base_type::state_type;
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

  using System::System;

  value_type velocity(state_type const& s) const { return s[0]; }
  state_type flux(state_type const& s) const {
    auto const& u = s[0];
    return state_type{0.5 * u * u};
  }
  state_type wave_speeds(state_type const& s) const { return s; }
  bool admissible(state_type const&) const { return true; }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }

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

/* ELASTIC COMPRESSION WAVES IN A SOLID */
struct LinearElasticPWaves : System<LinearElasticPWaves, fivo::state<double, 2>> {
  using base_type = System<LinearElasticPWaves, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
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

  LinearElasticPWaves(Mesh const& mesh,
                      std::shared_ptr<BC> const& left_bc,
                      std::shared_ptr<BC> const& right_bc,
                      value_type const& rho,
                      value_type const& lambda,
                      value_type const& mu)
    : base_type(mesh, left_bc, right_bc), m_rho(rho), m_lambda(lambda), m_mu(mu) {}

  value_type m_rho, m_lambda, m_mu;

  state_type flux(state_type const& s) const {
    return {- s[1], - (m_lambda + 2 * m_mu) / m_rho * s[0]};
  }
  state_type wave_speeds(state_type const& /* s */) const {
    auto const c = std::sqrt((m_lambda + 2 * m_mu) / m_rho);
    return state_type{-c, c};
  }
  bool admissible(state_type const&) const { return true; }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }
};

struct LinearElasticSWaves : System<LinearElasticSWaves, fivo::state<double, 2>> {
  using base_type = System<LinearElasticSWaves, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  LinearElasticSWaves(Mesh const& mesh,
                      std::shared_ptr<BC> const& left_bc,
                      std::shared_ptr<BC> const& right_bc,
                      value_type const& rho,
                      value_type const& mu)
    : base_type(mesh, left_bc, right_bc), m_rho(rho), m_mu(mu) {}

  value_type m_rho, m_mu;

  state_type flux(state_type const& s) const {
    return {- 0.5 * s[1], - 2 * m_mu / m_rho * s[0]};
  }
  state_type wave_speeds(state_type const& /* s */) const {
    auto const c = std::sqrt(m_mu / m_rho);
    return state_type{-c, c};
  }
  bool admissible(state_type const&) const { return true; }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }
};

/* LINEAR ACOUSTICS EQUATIONS : PRESSURE-VELOCITY FORMULATION */
struct LinearAcousticsPressure
  : System<LinearAcousticsPressure, fivo::state<double, 2>>,
    HasDensity<LinearAcousticsPressure, fivo::state<double, 2>>,
    HasPressure<LinearAcousticsPressure, fivo::state<double, 2>>,
    HasVelocity<LinearAcousticsPressure, fivo::state<double, 2>>,
    HasRiemannSolver<LinearAcousticsPressure, fivo::state<double, 2>> {

  using base_type = System<LinearAcousticsPressure, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
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

  LinearAcousticsPressure(Mesh const& mesh,
                          std::shared_ptr<BC> const& left_bc,
                          std::shared_ptr<BC> const& right_bc,
                          value_type const& c,
                          value_type const& r0,
                          value_type const& u0 = 0)
    : base_type(mesh, left_bc, right_bc), m_c(c), m_r0(r0), m_u0(u0) {}

  value_type m_c, m_r0, m_u0;

  value_type density(state_type const& s) const { return s[0] / (m_c * m_c); }
  value_type pressure(state_type const& s) const { return s[0]; }
  value_type velocity(state_type const& s) const { return s[1]; }
  state_type flux(state_type const& s) const {
    return {m_u0 * s[0] + m_r0 * m_c * m_c * s[1],
            s[0] / m_r0 + m_u0 * s[1]};
  }
  state_type wave_speeds(state_type const& /* s */) const {
    return state_type{m_u0 - m_c, m_u0 + m_c};
  }
  bool admissible(state_type const&) const { return true; }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }

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
};

/* LINEAR ACOUSTICS EQUATIONS : DENSITY-VELOCITY FORMULATION */
struct LinearAcousticsDensity
  : System<LinearAcousticsDensity, fivo::state<double, 2>>,
    HasDensity<LinearAcousticsDensity, fivo::state<double, 2>>,
    HasPressure<LinearAcousticsDensity, fivo::state<double, 2>>,
    HasVelocity<LinearAcousticsDensity, fivo::state<double, 2>>,
    HasRiemannSolver<LinearAcousticsDensity, fivo::state<double, 2>> {

  using base_type = System<LinearAcousticsDensity, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
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

  LinearAcousticsDensity(Mesh const& mesh,
                         std::shared_ptr<BC> const& left_bc,
                         std::shared_ptr<BC> const& right_bc,
                         value_type const& c,
                         value_type const& r0,
                         value_type const& u0 = 0)
    : base_type(mesh, left_bc, right_bc), m_c(c), m_r0(r0), m_u0(u0) {}

  value_type m_c, m_r0, m_u0;

  value_type density(state_type const& s) const { return s[0]; }
  value_type pressure(state_type const& s) const { return m_c * m_c * s[0]; }
  value_type velocity(state_type const& s) const { return s[1]; }
  state_type flux(state_type const& s) const {
    return {m_u0 * s[0] + m_r0 * s[1],
            m_c * m_c / m_r0 * s[0] + m_u0 * s[1]};
  }
  state_type wave_speeds(state_type const& /* s */) const {
    return state_type{m_u0 - m_c, m_u0 + m_c};
  }
  bool admissible(state_type const&) const { return true; }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& ul = left[1];
    auto const& rr = right[0];
    auto const& ur = right[1];
    auto const rstar = 0.5 * (rl + rr - (ur - ul) * m_r0 / m_c);
    auto const ustar = 0.5 * (rl + ur - (rr - rl) * m_c / m_r0);
    auto const wr = m_u0 - m_c;
    auto const wl = m_u0 + m_c;
    // Solution of the Riemann problem
    auto sol = [=] (value_type const& xt) {
                 if (xt < wl) return left;
                 if (xt > wr) return right;
                 return state_type{rstar, ustar};
               };
    return sol;
  }
};

/* SHALLOW WATER EQUATIONS */
struct SWE : System<SWE, fivo::state<double, 2>>,
             HasVelocity<SWE, fivo::state<double, 2>>,
             HasRiemannSolver<SWE, fivo::state<double, 2>> {
  using base_type = System<SWE, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
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

  SWE(Mesh const& mesh,
      std::shared_ptr<BC> const& left_bc,
      std::shared_ptr<BC> const& right_bc,
      value_type const& grav,
      std::shared_ptr<FrictionModel> const& fmodel = std::make_shared<NoFriction>(0))
    : base_type(mesh, left_bc, right_bc), m_grav(grav), m_fmodel(fmodel)
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

  value_type velocity(state_type const& s) const { return s[1] / s[0]; }
  state_type flux(state_type const& s) const {
    auto const& [h, j] = s;
    return state_type{j, j * j / h + 0.5 * m_grav * h * h};
  }
  state_type wave_speeds(state_type const& s) const {
    auto const h = s[0];
    auto const u = s[1] / h;
    auto const c = std::sqrt(m_grav * h);
    return state_type{u - c, u + c};
  }
  bool admissible(state_type const& s) const { return s[0] > 0; }
  state_type prim_to_cons(state_type const& s) const {
    auto const& [h, u] = s;
    return state_type{h, h * u};
  }
  state_type cons_to_prim(state_type const& s) const {
    auto const& [h, j] = s;
    return state_type{h, j / h};
  }

  global_state_type source(Mesh const& mesh, value_type const& t,
                           global_state_type const& X) const override {
    global_state_type S(mesh.nx());
    for (int i = 0; i < mesh.nx(); ++i) {
      S[i] = state_type{0., - m_gdz[i] * X[i][0]} + m_fmodel->apply(m_grav, X[i]);
    }
    return S;
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& hl = left[0];
    auto const& hr = right[0];
    auto const ul = velocity(left);
    auto const ur = velocity(right);
    auto const cl = std::sqrt(m_grav * hl);
    auto const cr = std::sqrt(m_grav * hr);

    // Find hstar and ustar in the star region
    auto const du = ur - ul;
    auto const fv =
      [=] (auto const h, auto const hk) { return std::sqrt(0.5 * m_grav * (1 / h + 1 / hk)); };

    auto const fk =
      [=] (auto const h, auto const hk) {
        if (h > hk) return (h - hk) * fv(h, hk);
        else return 2 * (std::sqrt(m_grav * h) - std::sqrt(m_grav * hk));
      };
    auto const dfk =
      [=] (auto const h, auto const hk) {
        if (h > hk) return fv(h, hk) + m_grav * (h - hk) * (hk - h) / (2 * h * h * fv(h, hk));
        else return std::sqrt(m_grav / h);
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
    auto const cstar = std::sqrt(m_grav * hstar);

    // Left shock speed
    auto const sl = ul - std::sqrt(0.5 * m_grav * hstar * hl * (hstar + hl)) / hl;

    // Left head and tail rarefaction speeds
    auto const shl = ul - cl;
    auto const stl = ustar - cstar;

    // Right shock speed
    auto const sr = ur + std::sqrt(0.5 * m_grav * hstar * hr * (hstar + hr)) / hr;

    // Right head and tail rarefaction speeds
    auto const shr = ur + cr;
    auto const str = ustar + cstar;

    // Solution in the star region\fan
    auto const lrstar =
      [=] (value_type const&) { return prim_to_cons(state_type{hstar, ustar}); };

    // Solution in the fan region for left rarefaction
    auto const lfan =
      [=] (value_type const& xt) {
        auto const hlfan = std::pow(2 * cl + ul - xt, 2) / (9 * m_grav);
        auto const ulfan = ul + 2. / 3. * (xt - ul + cl);
        return prim_to_cons(state_type{hlfan, ulfan});
      };

    // Solution in the fan region for right rarefaction
    auto const rfan =
      [=] (value_type const& xt) {
        auto const hrfan = std::pow(2 * cr - ur + xt, 2) / (9 * m_grav);
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
};

/* TELEGRAPH EQUATIONS */
struct Telegraph : System<Telegraph, fivo::state<double, 2>>,
                   HasRiemannSolver<Telegraph, fivo::state<double, 2>> {
  using base_type = System<Telegraph, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
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

  Telegraph(Mesh const& mesh,
            std::shared_ptr<BC> const& left_bc,
            std::shared_ptr<BC> const& right_bc,
            value_type const& velocity,
            value_type const& sigma)
    : base_type(mesh, left_bc, right_bc), m_velocity(velocity), m_sigma(sigma)
  {}

  value_type m_velocity;
  value_type m_sigma;

  auto velocity() const { return m_velocity; }
  auto sigma() const { return m_sigma; }

  state_type flux(state_type const& s) const {
    return state_type{m_velocity * s[0], -m_velocity * s[1]};
  }
  state_type wave_speeds(state_type const& /* s */) const {
    return state_type{-m_velocity, m_velocity};
  }
  bool admissible(state_type const& s) const { return (s[0] > 0 && s[1] > 0); }
  state_type prim_to_cons(state_type const& s) const { return s; }
  state_type cons_to_prim(state_type const& s) const { return s; }

  global_state_type source(Mesh const& mesh, value_type const& /* t */,
                           global_state_type const& X) const override {
    global_state_type src(mesh.nx());
    std::transform(X.begin(), X.end(), src.begin(),
                   [&] (auto const& s) {
                     return m_sigma * state_type{s[1] - s[0], s[0] - s[1]};
                   });
    return src;
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto exact = [&] (value_type const&) {
                   if (m_velocity > 0) return state_type{left[0], right[1]};
                   return state_type{left[1], right[0]};
                 };
    return exact;
  }
};

/* ISENTROPIC EULER EQUATIONS */
struct IsentropicEuler : System<IsentropicEuler, fivo::state<double, 2>>,
                         HasDensity<IsentropicEuler, fivo::state<double, 2>>,
                         HasVelocity<IsentropicEuler, fivo::state<double, 2>>,
                         HasPressure<IsentropicEuler, fivo::state<double, 2>> {
  using base_type = System<IsentropicEuler, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
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

  IsentropicEuler(Mesh const& mesh,
                  std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc,
                  value_type const& kappa,
                  value_type const& gamma)
    : base_type(mesh, left_bc, right_bc), m_kappa(kappa), m_gamma(gamma) {}

  value_type m_kappa, m_gamma;

  auto kappa() const { return m_kappa; }
  auto gamma() const { return m_gamma; }

  value_type density(state_type const& s) const { return s[0]; }
  value_type velocity(state_type const& s) const { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const {
    return m_kappa * std::pow(s[0], m_gamma);
  }
  state_type flux(state_type const& s) const {
    auto const p = pressure(s);
    auto const& [r, j] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p};
  }
  state_type wave_speeds(state_type const& s) const {
    auto const p = pressure(s);
    auto const& [r, u] = cons_to_prim(s);
    auto const c = std::sqrt(m_gamma * p / r);
    return state_type{u - c, u + c};
  }
  bool admissible(state_type const& s) const { return (s[0] > 0); }
  state_type prim_to_cons(state_type const& s) const {
    auto const& [r, u] = s;
    return state_type{r, r * u};
  }
  state_type cons_to_prim(state_type const& s) const {
    auto const& [r, j] = s;
    return state_type{r, j / r};
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& rr = right[0];
    auto const ul = velocity(left);
    auto const ur = velocity(right);
    auto const cl = std::sqrt(m_kappa * m_gamma * std::pow(rl, m_gamma - 1));
    auto const cr = std::sqrt(m_kappa * m_gamma * std::pow(rr, m_gamma - 1));

    // Find hstar and ustar in the star region
    auto const du = ur - ul;
    auto const fk =
      [=] (auto const r, auto const rk) {
        if (r > rk) {
          auto const p = m_kappa * std::pow(r, m_gamma - 1);
          auto const pk = m_kappa * std::pow(rk, m_gamma - 1);
          return (r - rk) * std::sqrt((pk - p) / (r * rk * (rk - r))); 
        }
        else {
          auto const cstar = m_kappa * m_gamma * std::pow(r, m_gamma - 1);
          return 2 / (m_gamma - 1) * (cstar - cr);
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
    auto const cstar = std::sqrt(m_kappa * m_gamma * std::pow(rstar, m_gamma - 1));

    // Left shock speed
    auto const sl = ul - 2 / (m_gamma - 1) * (cstar - cl);

    // Left head and tail rarefaction speeds
    auto const shl = ul - cl;
    auto const stl = ustar - cstar;

    // Right shock speed
    auto const sr = ur - 2 / (m_gamma - 1) * (cstar - cr);

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
};

/* ISOTHERMAL EULER EQUATIONS */
struct IsothermalEuler : System<IsothermalEuler, fivo::state<double, 2>>,
                         HasDensity<IsothermalEuler, fivo::state<double, 2>>,
                         HasVelocity<IsothermalEuler, fivo::state<double, 2>>,
                         HasPressure<IsothermalEuler, fivo::state<double, 2>> {
  using base_type = System<IsothermalEuler, fivo::state<double, 2>>;
  using state_type = typename base_type::state_type;
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

  IsothermalEuler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                  std::shared_ptr<BC> const& right_bc, value_type const& c)
    : base_type(mesh, left_bc, right_bc), m_c(c) {}

  value_type m_c;

  auto c() const { return m_c; }

  value_type density(state_type const& s) const { return s[0]; }
  value_type velocity(state_type const& s) const { return s[1] / s[0]; }
  value_type pressure(state_type const& s) const { return m_c * m_c * s[0]; }
  state_type flux(state_type const& s) const {
    auto const p = pressure(s);
    auto const& [r, j] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p};
  }
  state_type wave_speeds(state_type const& s) const {
    auto const u = s[1] / s[0];
    return state_type{u - m_c, u + m_c};
  }
  bool admissible(state_type const& s) const { return (s[0] > 0); }
  state_type prim_to_cons(state_type const& s) const {
    auto const& [r, u] = s;
    return state_type{r, r * u};
  }
  state_type cons_to_prim(state_type const& s) const {
    auto const& [r, j] = s;
    return state_type{r, j / r};
  }
};

/* GENERAL EULER EQUATIONS (ABSTRACT) */
template<typename Derived>
struct Euler : System<Derived, fivo::state<double, 3>>,
               HasDensity<Derived, fivo::state<double, 3>>,
               HasVelocity<Derived, fivo::state<double, 3>>,
               HasPressure<Derived, fivo::state<double, 3>> {
  using base_type = System<Derived, fivo::state<double, 3>>;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  /* BOUNDARY CONDITIONS */
  using BC = typename base_type::BC;
  struct BCReflective : BC {
    static inline auto make() { return std::make_shared<BCReflective>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      return state_type{in[0], -in[1], in[2]};
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

  Euler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
        std::shared_ptr<BC> const& right_bc, value_type const& gamma)
    : base_type(mesh, left_bc, right_bc), m_gamma(gamma) {}

  value_type m_gamma;

  auto gamma() const { return m_gamma; }

  value_type density(state_type const& s) const { return s[0]; }
  value_type velocity(state_type const& s) const { return s[1] / s[0]; }
  state_type flux(state_type const& s) const {
    auto const p = static_cast<Derived const*>(this)->pressure(s);
    auto const& [r, j, e] = s;
    auto const u = j / r;
    return state_type{j, r * u * u + p, (e + p) * u};
  }
  state_type wave_speeds(state_type const& s) const {
    auto const p = static_cast<Derived const*>(this)->pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_gamma * p / r);
    return state_type{u - c, u, u + c};
  }
  bool admissible(state_type const& s) const {
    auto const r = s[0];
    auto const e = s[2];
    return (r > 0 && e > 0);
  }
};

/* GENERAL EULER EQUATIONS WITH AN IDEAL GAS EOS */
struct IdealGasEuler : Euler<IdealGasEuler>,
                       HasRiemannSolver<IdealGasEuler,
                                        typename Euler<IdealGasEuler>::state_type> {
  using base_type = Euler;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  using Euler::Euler;

  value_type pressure(state_type const& s) const {
    auto const& [r, j, e] = s;
    return (m_gamma - 1) * (e - 0.5 * j * j / r);
  }
  value_type energy(value_type const& r, value_type const& u, value_type const& p) const {
    return p / (m_gamma - 1) + 0.5 * r * u * u;
  }
  state_type prim_to_cons(state_type const& s) const {
    auto const& [r, u, p] = s;
    return state_type{r, r * u, energy(r, u, p)};
  }
  state_type cons_to_prim(state_type const& s) const {
    auto const& [r, j, e] = s;
    return state_type{r, j / r, pressure(s)};
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& rr = right[0];
    auto const pl = pressure(left);
    auto const pr = pressure(right);
    auto const ul = velocity(left);
    auto const ur = velocity(right);
    auto const cl = std::sqrt(m_gamma * pl / rl);
    auto const cr = std::sqrt(m_gamma * pr / rr);

    // Find pstar and ustar in the star region
    auto const Al = 2. / (rl * (m_gamma + 1));
    auto const Ar = 2. / (rr * (m_gamma + 1));
    auto const Bl = pl * (m_gamma - 1) / (m_gamma + 1);
    auto const Br = pr * (m_gamma - 1) / (m_gamma + 1);
    auto const du = ur - ul;

    auto const fk =
      [=] (auto const p, auto const pk, auto const ck, auto const Ak, auto const Bk) {
        if (p > pk) return (p - pk) * std::sqrt(Ak / (p + Bk));
        else return 2 * ck / (m_gamma - 1) * (std::pow(p / pk, (m_gamma - 1) / (2 * m_gamma)) - 1);
      };
    auto const dfk =
      [=] (auto const p, auto const pk, auto const rk,
          auto const ck, auto const Ak, auto const Bk) {
        if (p > pk) return std::sqrt(Ak / (Bk + p)) * (1 - (p - pk) / (2 * (Bk + p)));
        else return std::pow(p / pk, - (m_gamma + 1) / (2 * m_gamma)) / (rk * ck);
      };
    auto const fl = [&] (auto const p) { return fk(p, pl, cl, Al, Bl); };
    auto const fr = [&] (auto const p) { return fk(p, pr, cr, Ar, Br); };
    auto const dfl = [&] (auto const p) { return dfk(p, pl, rl, cl, Al, Bl); };
    auto const dfr = [&] (auto const p) { return dfk(p, pr, rr, cr, Ar, Br); };
    auto const f = [&] (auto const p) { return fl(p) + fr(p) + du; };
    auto const df = [&] (auto const p) { return dfl(p) + dfr(p); };

    auto const tol = 1e-10;
    auto const ppv = 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rl + rr) * (cl + cr);
    auto const p0 = std::max(tol, ppv);

    auto const pstar = math::newton_raphson(p0, f, df, tol);
    auto const ustar = 0.5 * ((ul + ur) + fr(pstar) - fl(pstar));

    /* FIND rlstar AND rrstar (density in the star region) */
    auto const fac = (m_gamma - 1) / (m_gamma + 1);

    // Solution left of the contact
    auto const rlstar = (pstar > pl)
      ? rl * (pstar / pl + fac) / ((pstar / pl) * fac + 1) // left shock
      : rl * std::pow(pstar / pl, 1 / m_gamma);            // left rarefaction

    // Solution right of the contact
    auto const rrstar = (pstar > pr)
      ? rr * (pstar / pr + fac) / ((pstar / pr) * fac + 1) // right shock
      : rr * std::pow(pstar / pr, 1 / m_gamma);            // right rarefaction

    // Left shock speed
    auto const sl = ul - cl * std::sqrt(((m_gamma + 1) * pstar / pl + (m_gamma - 1)) / (2 * m_gamma));
    // Left head and tail rarefaction speeds
    auto const clstar = cl * std::pow(pstar / pl, (m_gamma - 1) / (2 * m_gamma));
    auto const shl = ul - cl;
    auto const stl = ustar - clstar;

    // Right shock speed
    auto const sr = ur + cr * std::sqrt(((m_gamma + 1) * pstar / pr + (m_gamma - 1)) / (2 * m_gamma));
    // Right head and tail rarefaction speeds
    auto const crstar = cr * std::pow(pstar / pr, (m_gamma - 1) / (2 * m_gamma));
    auto const shr = ur + cr;
    auto const str = ustar + crstar;

    // Solution in the left star region\fan
    auto const lstar =
      [=] (value_type const&) { return prim_to_cons(state_type{rlstar, ustar, pstar}); };

    // Solution in the right star region\fan
    auto const rstar =
      [=] (value_type const&) { return prim_to_cons(state_type{rrstar, ustar, pstar}); };

    // Solution in the fan region for left rarefaction
    auto const lfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 + (m_gamma - 1) * (ul - xt) / cl) / (m_gamma + 1);
        auto const rlfan = rl * std::pow(tmp, 2 / (m_gamma - 1));
        auto const ulfan = 2 / (m_gamma + 1) * (cl + (m_gamma - 1) * ul / 2 + xt);
        auto const plfan = pl * std::pow(tmp, 2 * m_gamma / (m_gamma - 1));
        return prim_to_cons(state_type{rlfan, ulfan, plfan});
      };

    // Solution in the fan region for right rarefaction
    auto const rfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 - (m_gamma - 1) * (ur - xt) / cr) / (m_gamma + 1);
        auto const rrfan = rr * std::pow(tmp, 2 / (m_gamma - 1));
        auto const urfan = 2 / (m_gamma + 1) * (-cr + (m_gamma - 1) * ur / 2 + xt);
        auto const prfan = pr * std::pow(tmp, 2 * m_gamma / (m_gamma - 1));
        return prim_to_cons(state_type{rrfan, urfan, prfan});
      };

    auto const sol =
      [=] (value_type const& xt) {
        if (xt < ustar) { // Left of contact
          if (pstar > pl) { return (xt < sl) ? left : lstar(xt); }
          else { return (xt < shl) ? left : (xt > stl) ? lstar(xt) : lfan(xt); }
        }
        else { // Right of contact
          if (pstar > pr) { return (xt > sr) ? right : rstar(xt);}
          else { return (xt > shr) ? right : (xt < str) ? rstar(xt) : rfan(xt); }
        }
      };
    return sol;
  }
};

/* GENERAL EULER EQUATIONS WITH A COVOLUME GAS EOS */
struct CovolumeGasEuler : Euler<CovolumeGasEuler>,
                          HasRiemannSolver<CovolumeGasEuler,
                                           typename Euler<CovolumeGasEuler>::state_type> {
  using base_type = Euler;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  CovolumeGasEuler(Mesh const& mesh,
                   std::shared_ptr<BC> const& left_bc,
                   std::shared_ptr<BC> const& right_bc,
                   value_type const& gamma, value_type const& b)
    : base_type(mesh, left_bc, right_bc, gamma), m_b(b) {}

  value_type m_b;

  value_type pressure(state_type const& s) const {
    auto const& [r, j, e] = s;
    return (m_gamma - 1) * (e - 0.5 * j * j / r) / (1 - m_b * r);
  }
  value_type energy(value_type const& r, value_type const& u, value_type const& p) const {
    return p * (1 - m_b * r) / (m_gamma - 1) + 0.5 * r * u * u;
  }
  state_type prim_to_cons(state_type const& s) const {
    auto const& [r, u, p] = s;
    return state_type{r, r * u, energy(r, u, p)};
  }
  state_type cons_to_prim(state_type const& s) const {
    auto const& [r, j, e] = s;
    return state_type{r, j / r, pressure(s)};
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& rl = left[0];
    auto const& rr = right[0];
    auto const pl = pressure(left);
    auto const pr = pressure(right);
    auto const ul = velocity(left);
    auto const ur = velocity(right);
    auto const cl = std::sqrt(m_gamma * pl / (rl * (1 - m_b * rl)));
    auto const cr = std::sqrt(m_gamma * pr / (rr * (1 - m_b * rr)));

    // Find pstar and ustar in the star region
    auto const Al = 2. * (1 - m_b * rl) / (rl * (m_gamma + 1));
    auto const Ar = 2. * (1 - m_b * rr) / (rr * (m_gamma + 1));
    auto const Bl = pl * (m_gamma - 1) / (m_gamma + 1);
    auto const Br = pr * (m_gamma - 1) / (m_gamma + 1);
    auto const du = ur - ul;

    auto const fk =
      [=] (auto const p, auto const pk, auto const ck, auto const Ak, auto const Bk) {
        if (p > pk) return (p - pk) * std::sqrt(Ak / (p + Bk));
        else return 2 * ck * (1 - m_b * pk) / (m_gamma - 1) * (std::pow(p / pk, (m_gamma - 1) / (2 * m_gamma)) - 1);
      };
    auto const dfk =
      [=] (auto const p, auto const pk, auto const rk,
          auto const ck, auto const Ak, auto const Bk) {
        if (p > pk) return std::sqrt(Ak / (Bk + p)) * (1 - (p - pk) / (2 * (Bk + p)));
        else return (1 - m_b * pk) * std::pow(p / pk, - (m_gamma + 1) / (2 * m_gamma)) / (rk * ck);
      };
    auto const fl = [&] (auto const p) { return fk(p, pl, cl, Al, Bl); };
    auto const fr = [&] (auto const p) { return fk(p, pr, cr, Ar, Br); };
    auto const dfl = [&] (auto const p) { return dfk(p, pl, rl, cl, Al, Bl); };
    auto const dfr = [&] (auto const p) { return dfk(p, pr, rr, cr, Ar, Br); };
    auto const f = [&] (auto const p) { return fl(p) + fr(p) + du; };
    auto const df = [&] (auto const p) { return dfl(p) + dfr(p); };

    auto const tol = 1e-10;
    auto const ppv = 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rl + rr) * (cl + cr);
    auto const p0 = std::max(tol, ppv);

    auto const pstar = math::newton_raphson(p0, f, df, tol);
    auto const ustar = 0.5 * ((ul + ur) + fr(pstar) - fl(pstar));

    /* FIND rlstar AND rrstar (density in the star region) */
    auto const fac = (m_gamma - 1) / (m_gamma + 1);

    // Solution left of the contact
    auto const rlstar = (pstar > pl)
      ? rl * (pstar / pl + fac) / ((pstar / pl) * fac + 1) // left shock
      : rl * std::pow(pstar / pl, 1 / m_gamma);            // left rarefaction

    // Solution right of the contact
    auto const rrstar = (pstar > pr)
      ? rr * (pstar / pr + fac) / ((pstar / pr) * fac + 1) // right shock
      : rr * std::pow(pstar / pr, 1 / m_gamma);            // right rarefaction

    // Left shock speed
    auto const sl = ul - cl * std::sqrt(((m_gamma + 1) * pstar / pl + (m_gamma - 1)) / (2 * m_gamma));
    // Left head and tail rarefaction speeds
    auto const clstar = cl * std::pow(pstar / pl, (m_gamma - 1) / (2 * m_gamma));
    auto const shl = ul - cl;
    auto const stl = ustar - clstar;

    // Right shock speed
    auto const sr = ur + cr * std::sqrt(((m_gamma + 1) * pstar / pr + (m_gamma - 1)) / (2 * m_gamma));
    // Right head and tail rarefaction speeds
    auto const crstar = cr * std::pow(pstar / pr, (m_gamma - 1) / (2 * m_gamma));
    auto const shr = ur + cr;
    auto const str = ustar + crstar;

    // Solution in the left star region\fan
    auto const lstar =
      [=] (value_type const&) { return prim_to_cons(state_type{rlstar, ustar, pstar}); };

    // Solution in the right star region\fan
    auto const rstar =
      [=] (value_type const&) { return prim_to_cons(state_type{rrstar, ustar, pstar}); };

    // Solution in the fan region for left rarefaction
    auto const lfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 + (m_gamma - 1) * (ul - xt) / cl) / (m_gamma + 1);
        auto const rlfan = rl * std::pow(tmp, 2 / (m_gamma - 1));
        auto const ulfan = 2 / (m_gamma + 1) * (cl + (m_gamma - 1) * ul / 2 + xt);
        auto const plfan = pl * std::pow(tmp, 2 * m_gamma / (m_gamma - 1));
        return prim_to_cons(state_type{rlfan, ulfan, plfan});
      };

    // Solution in the fan region for right rarefaction
    auto const rfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 - (m_gamma - 1) * (ur - xt) / cr) / (m_gamma + 1);
        auto const rrfan = rr * std::pow(tmp, 2 / (m_gamma - 1));
        auto const urfan = 2 / (m_gamma + 1) * (-cr + (m_gamma - 1) * ur / 2 + xt);
        auto const prfan = pr * std::pow(tmp, 2 * m_gamma / (m_gamma - 1));
        return prim_to_cons(state_type{rrfan, urfan, prfan});
      };

    auto const sol =
      [=] (value_type const& xt) {
        if (xt < ustar) { // Left of contact
          if (pstar > pl) { return (xt < sl) ? left : lstar(xt); }
          else { return (xt < shl) ? left : (xt > stl) ? lstar(xt) : lfan(xt); }
        }
        else { // Right of contact
          if (pstar > pr) { return (xt > sr) ? right : rstar(xt);}
          else { return (xt > shr) ? right : (xt < str) ? rstar(xt) : rfan(xt); }
        }
      };
    return sol;
  }
};

/* GENERAL EULER EQUATIONS WITH A STIFFENED GAS EOS */
struct StiffenedGasEuler : Euler<StiffenedGasEuler> {
  using base_type = Euler;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  StiffenedGasEuler(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
                    std::shared_ptr<BC> const& right_bc, value_type const& gamma,
                    value_type const& p0)
    : Euler(mesh, left_bc, right_bc, gamma), m_p0(p0) {}

  value_type m_p0;
  auto p0() const { return m_p0; }

  value_type pressure(state_type const& s) const {
    auto const& [r, j, e] = s;
    return (m_gamma - 1) * (e - 0.5 * j * j / r) - m_gamma * m_p0;
  }
  value_type energy(value_type const& r, value_type const& u, value_type const& p) const {
    return (p + m_gamma * m_p0) / (m_gamma - 1) + 0.5 * r * u * u;
  }
  state_type prim_to_cons(state_type const& s) const {
    auto const& [r, u, p] = s;
    return state_type{r, r * u, energy(r, u, p)};
  }
  state_type cons_to_prim(state_type const& s) const {
    auto const& [r, j, e] = s;
    return state_type{r, j / r, pressure(s)};
  }
};

/* GENERAL EULER EQUATIONS WITH PASSIVE SCALARS (ABSTRACT) */
template<typename Derived, std::size_t NumPassiveScalars = 0>
struct EulerP
  : System<EulerP<Derived, NumPassiveScalars>, fivo::state<double, 3 + NumPassiveScalars>>,
  HasDensity<EulerP<Derived, NumPassiveScalars>, fivo::state<double, 3 + NumPassiveScalars>>,
  HasVelocity<EulerP<Derived, NumPassiveScalars>, fivo::state<double, 3 + NumPassiveScalars>>,
  HasPressure<EulerP<Derived, NumPassiveScalars>, fivo::state<double, 3 + NumPassiveScalars>> {

  using base_type = System<EulerP<Derived, NumPassiveScalars>,
                           fivo::state<double, 3 + NumPassiveScalars>>;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  static constexpr std::size_t num_passive_scalars = NumPassiveScalars;

  /* BOUNDARY CONDITIONS */
  using BC = typename base_type::BC;
  struct BCReflective : BC {
    static inline auto make() { return std::make_shared<BCReflective>(); }
    state_type compute(Mesh const& mesh, value_type const& t,
                       state_type const& in, state_type const& in_opbound) const final {
      state_type ghost = {in[0], -in[1], in[2]};
      for (std::size_t i = 0; i < num_passive_scalars; ++i) { ghost[3 + i] = -in[3 + i]; }
      return ghost;
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

  EulerP(Mesh const& mesh, std::shared_ptr<BC> const& left_bc,
         std::shared_ptr<BC> const& right_bc, value_type const& gamma)
    : base_type(mesh, left_bc, right_bc), m_gamma(gamma) {}

  value_type m_gamma;

  auto gamma() const { return m_gamma; }

  value_type pressure(state_type const& s) const {
    return static_cast<Derived const*>(this)->pressure(s);
  }
  value_type density(state_type const& s) const { return s[0]; }
  value_type velocity(state_type const& s) const { return s[1] / s[0]; }
  state_type flux(state_type const& s) const {
    auto const p = this->pressure(s);
    auto const& r = s[0];
    auto const& j = s[1];
    auto const& e = s[2];
    auto const u = j / r;
    state_type fl = {j, r * u * u + p, (e + p) * u};
    for (std::size_t i = 0; i < num_passive_scalars; ++i) { fl[3 + i] = u * s[3 + i]; }
    return fl;
  }
  state_type wave_speeds(state_type const& s) const {
    auto const p = this->pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_gamma * p / r);
    state_type ws = {u - c, u, u + c};
    for (std::size_t i = 0; i < num_passive_scalars; ++i) { ws[3 + i] = u; }
    return ws;
  }
  bool admissible(state_type const& s) const {
    auto const r = s[0];
    auto const e = s[2];
    return (r > 0 && e > 0);
  }
};

/* GENERAL EULER EQUATIONS WITH AN IDEAL GAS EOS */
template<std::size_t NumPassiveScalars = 0>
struct IdealGasEulerP : EulerP<IdealGasEulerP<NumPassiveScalars>, NumPassiveScalars>,
  HasRiemannSolver<IdealGasEulerP<NumPassiveScalars>,
                   typename EulerP<IdealGasEulerP<NumPassiveScalars>,
                                   NumPassiveScalars>::state_type> {

  using base_type = EulerP<IdealGasEulerP<NumPassiveScalars>, NumPassiveScalars>;
  using state_type = typename base_type::state_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  using base_type::EulerP;
  using base_type::num_passive_scalars;

  value_type pressure(state_type const& s) const {
    auto const& r = s[0];
    auto const& j = s[1];
    auto const& e = s[2];
    return (this->m_gamma - 1) * (e - 0.5 * j * j / r);
  }
  value_type energy(value_type const& r, value_type const& u, value_type const& p) const {
    return p / (this->m_gamma - 1) + 0.5 * r * u * u;
  }
  state_type prim_to_cons(state_type const& s) const {
    auto const& r = s[0];
    auto const& u = s[1];
    auto const& p = s[2];
    state_type cons = {r, r * u, energy(r, u, p)};
    for (std::size_t i = 0; i < num_passive_scalars; ++i) { cons[3 + i] = r * s[3 + i]; }
    return cons;
  }
  state_type cons_to_prim(state_type const& s) const {
    auto const& r = s[0];
    auto const& j = s[1];
    state_type prim = {r, j / r, pressure(s)};
    for (std::size_t i = 0; i < num_passive_scalars; ++i) { prim[3 + i] = s[3 + i] / r; }
    return prim;
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& gamma = this->m_gamma;
    auto const pleft = cons_to_prim(left);
    auto const pright = cons_to_prim(right);
    auto const& rl = pleft[0];
    auto const& ul = pleft[1];
    auto const& pl = pleft[2];
    auto const& rr = pright[0];
    auto const& ur = pright[1];
    auto const& pr = pright[2];
    auto const cl = std::sqrt(gamma * pl / rl);
    auto const cr = std::sqrt(gamma * pr / rr);

    // Find pstar and ustar in the star region
    auto const Al = 2. / (rl * (gamma + 1));
    auto const Ar = 2. / (rr * (gamma + 1));
    auto const Bl = pl * (gamma - 1) / (gamma + 1);
    auto const Br = pr * (gamma - 1) / (gamma + 1);
    auto const du = ur - ul;

    auto const fk =
      [=] (auto const p, auto const pk, auto const ck, auto const Ak, auto const Bk) {
        if (p > pk) return (p - pk) * std::sqrt(Ak / (p + Bk));
        else return 2 * ck / (gamma - 1) * (std::pow(p / pk, (gamma - 1) / (2 * gamma)) - 1);
      };
    auto const dfk =
      [=] (auto const p, auto const pk, auto const rk,
          auto const ck, auto const Ak, auto const Bk) {
        if (p > pk) return std::sqrt(Ak / (Bk + p)) * (1 - (p - pk) / (2 * (Bk + p)));
        else return std::pow(p / pk, - (gamma + 1) / (2 * gamma)) / (rk * ck);
      };
    auto const fl = [&] (auto const p) { return fk(p, pl, cl, Al, Bl); };
    auto const fr = [&] (auto const p) { return fk(p, pr, cr, Ar, Br); };
    auto const dfl = [&] (auto const p) { return dfk(p, pl, rl, cl, Al, Bl); };
    auto const dfr = [&] (auto const p) { return dfk(p, pr, rr, cr, Ar, Br); };
    auto const f = [&] (auto const p) { return fl(p) + fr(p) + du; };
    auto const df = [&] (auto const p) { return dfl(p) + dfr(p); };

    auto const tol = 1e-10;
    auto const ppv = 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rl + rr) * (cl + cr);
    auto const p0 = std::max(tol, ppv);

    auto const pstar = math::newton_raphson(p0, f, df, tol);
    auto const ustar = 0.5 * ((ul + ur) + fr(pstar) - fl(pstar));

    /* FIND rlstar AND rrstar (density in the star region) */
    auto const fac = (gamma - 1) / (gamma + 1);

    // Solution left of the contact
    auto const rlstar = (pstar > pl)
      ? rl * (pstar / pl + fac) / ((pstar / pl) * fac + 1) // left shock
      : rl * std::pow(pstar / pl, 1 / gamma);            // left rarefaction

    // Solution right of the contact
    auto const rrstar = (pstar > pr)
      ? rr * (pstar / pr + fac) / ((pstar / pr) * fac + 1) // right shock
      : rr * std::pow(pstar / pr, 1 / gamma);            // right rarefaction

    // Left shock speed
    auto const sl = ul - cl * std::sqrt(((gamma + 1) * pstar / pl + (gamma - 1)) / (2 * gamma));
    // Left head and tail rarefaction speeds
    auto const clstar = cl * std::pow(pstar / pl, (gamma - 1) / (2 * gamma));
    auto const shl = ul - cl;
    auto const stl = ustar - clstar;

    // Right shock speed
    auto const sr = ur + cr * std::sqrt(((gamma + 1) * pstar / pr + (gamma - 1)) / (2 * gamma));
    // Right head and tail rarefaction speeds
    auto const crstar = cr * std::pow(pstar / pr, (gamma - 1) / (2 * gamma));
    auto const shr = ur + cr;
    auto const str = ustar + crstar;

    // Solution in the left star region\fan
    auto const lstar =
      [=] (value_type const&) {
        state_type prim = {rlstar, ustar, pstar};
        for (std::size_t i = 0; i < num_passive_scalars; ++i) { prim[3 + i] = pleft[3 + i]; }
        return prim_to_cons(prim);
      };

    // Solution in the right star region\fan
    auto const rstar =
      [=] (value_type const&) {
        state_type prim = {rrstar, ustar, pstar};
        for (std::size_t i = 0; i < num_passive_scalars; ++i) { prim[3 + i] = pright[3 + i]; }
        return prim_to_cons(prim);
      };

    // Solution in the fan region for left rarefaction
    auto const lfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 + (gamma - 1) * (ul - xt) / cl) / (gamma + 1);
        auto const rlfan = rl * std::pow(tmp, 2 / (gamma - 1));
        auto const ulfan = 2 / (gamma + 1) * (cl + (gamma - 1) * ul / 2 + xt);
        auto const plfan = pl * std::pow(tmp, 2 * gamma / (gamma - 1));
        state_type prim = {rlfan, ulfan, plfan};
        for (std::size_t i = 0; i < num_passive_scalars; ++i) { prim[3 + i] = pleft[3 + i]; }
        return prim_to_cons(prim);
      };

    // Solution in the fan region for right rarefaction
    auto const rfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 - (gamma - 1) * (ur - xt) / cr) / (gamma + 1);
        auto const rrfan = rr * std::pow(tmp, 2 / (gamma - 1));
        auto const urfan = 2 / (gamma + 1) * (-cr + (gamma - 1) * ur / 2 + xt);
        auto const prfan = pr * std::pow(tmp, 2 * gamma / (gamma - 1));
        state_type prim = {rrfan, urfan, prfan};
        for (std::size_t i = 0; i < num_passive_scalars; ++i) { prim[3 + i] = pright[3 + i]; }
        return prim_to_cons(prim);
      };

    auto const sol =
      [=] (value_type const& xt) {
        if (xt < ustar) { // Left of contact
          if (pstar > pl) { return (xt < sl) ? left : lstar(xt); }
          else { return (xt < shl) ? left : (xt > stl) ? lstar(xt) : lfan(xt); }
        }
        else { // Right of contact
          if (pstar > pr) { return (xt > sr) ? right : rstar(xt);}
          else { return (xt > shr) ? right : (xt < str) ? rstar(xt) : rfan(xt); }
        }
      };
    return sol;
  }
};

} // namespace system

} // namespace fivo

#endif // FIVO_SYSTEM_H
