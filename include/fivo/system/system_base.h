#ifndef FIVO_SYSTEM_BASE_H
#define FIVO_SYSTEM_BASE_H

#include <fivo/state.h>
#include <fivo/mesh.h>
#include <fivo/vector.h>
#include <fivo/math.h>
#include <fivo/traits.h>

#include <memory>

namespace fivo {

namespace system {

template<typename State>
struct HasVelocity {
  using state_type = State;
  using global_state_type = global_state<state_type>;
  using value_type = typename state_type::value_type;

  virtual value_type velocity(state_type const& in) const = 0;
  auto velocities(global_state_type const& in) const {
    vector<value_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (auto const& s) { return velocity(s); });
    return out;
  }
};

template<typename State>
struct HasDensity {
  using state_type = State;
  using global_state_type = global_state<state_type>;
  using value_type = typename state_type::value_type;

  virtual value_type density(state_type const& in) const = 0;
  auto densities(global_state_type const& in) const {
    vector<value_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (auto const& s) { return density(s); });
    return out;
  }
};

template<typename State>
struct HasPressure {
  using state_type = State;
  using global_state_type = global_state<state_type>;
  using value_type = typename state_type::value_type;

  virtual value_type pressure(state_type const& in) const = 0;
  auto pressures(global_state_type const& in) const {
    vector<value_type> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(),
                   [&] (auto const& s) { return pressure(s); });
    return out;
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
    traits::is_derived_v<System, system::HasVelocity<typename System::state_type>>;
};
template<typename System>
constexpr bool has_velocity_v = has_velocity<System>::value;

template<typename System>
struct has_density {
  static constexpr bool value =
    traits::is_derived_v<System, system::HasDensity<typename System::state_type>>;
};
template<typename System>
constexpr bool has_density_v = has_density<System>::value;

template<typename System>
struct has_pressure {
  static constexpr bool value =
    traits::is_derived_v<System, system::HasPressure<typename System::state_type>>;
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

template<typename State>
struct SystemBase {
  using state_type = State;
  using flux_type = fivo::flux<State, 1>;
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

  SystemBase(Mesh const& mesh,
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
  virtual flux_type flux(state_type const& in) const = 0;

  /**
   * \brief Computes the wave speeds (i.e. the eigenvalues of the flux jacobian)
   * for a given local state.
   *
   * \param[in] in The input local state
   * \return All wave speeds
   */
  virtual state_type wave_speeds(state_type const& in) const = 0;

  /**
   * \brief Determines if a local state is admissible
   *
   * \param[in] in The input local state
   * \return True or False depending on the admissibility conditions
   */
  virtual bool admissible(state_type const& in) const = 0;

  /**
   * \brief Transforms a local state in conservative variables to primitives variables.
   *
   * For some systems, this may be the identity function.
   *
   * \param[in] in The input local state in conservative variables
   * \return The same physical local state in primitive variables formulation
   */
  virtual state_type cons_to_prim(state_type const& in) const = 0;

  /**
   * \brief Transforms a local state in primitive variables to conservative variables.
   *
   * For some systems, this may be the identity function.
   *
   * \param[in] in The input local state in primitive variables
   * \return The same physical local state in conservative variables formulation
   */
  virtual state_type prim_to_cons(state_type const& in) const = 0;

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

} // namespace system

} // namespace fivo

#endif // FIVO_SYSTEM_BASE_H
