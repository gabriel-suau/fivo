#ifndef FIVO_EULER_BASE_H
#define FIVO_EULER_BASE_H

#include <fivo/system/system_base.h>

namespace fivo {

namespace system {

// class EulerBase : public SystemBase<fivo::state<double, 3>>,
//                   public HasDensity<fivo::state<double, 3>>,
//                   public HasVelocity<fivo::state<double, 3>>,
//                   public HasPressure<fivo::state<double, 3>> {
// public:
//   using base_type = SystemBase<fivo::state<double, 3>>;
//   using state_type = typename base_type::state_type;
//   using flux_type = typename base_type::flux_type;
//   using global_state_type = typename base_type::global_state_type;
//   using value_type = typename base_type::value_type;

//   /* BOUNDARY CONDITIONS */
//   using BC = typename base_type::BC;
//   struct BCReflective : BC {
//     static inline auto make() { return std::make_shared<BCReflective>(); }
//     state_type compute(Mesh const& mesh, value_type const& t,
//                        state_type const& in, state_type const& in_opbound) const final {
//       return state_type{in[0], -in[1], in[2]};
//     }
//   };
//   struct BCTransmissive : BC {
//     static inline auto make() { return std::make_shared<BCTransmissive>(); }
//     state_type compute(Mesh const& mesh, value_type const& t,
//                        state_type const& in, state_type const& in_opbound) const final {
//       return in;
//     }
//   };
//   struct BCPeriodic : BC {
//     static inline auto make() { return std::make_shared<BCPeriodic>(); }
//     state_type compute(Mesh const& mesh, value_type const& t,
//                        state_type const& in, state_type const& in_opbound) const final {
//       return in_opbound;
//     }
//   };

//   struct Params { value_type gamma; };

//   EulerBase(Mesh const& mesh,
//             std::shared_ptr<BC> const& left_bc,
//             std::shared_ptr<BC> const& right_bc,
//             Params const& params)
//     : base_type(mesh, left_bc, right_bc), m_params(params) {}

//   value_type density(state_type const& s) const override { return s[0]; }
//   value_type velocity(state_type const& s) const override { return s[1] / s[0]; }
//   flux_type flux(state_type const& s) const override {
//     auto const p = pressure(s);
//     auto const& [r, j, e] = s;
//     auto const u = j / r;
//     return flux_type{{{{j}, {r * u * u + p}, {(e + p) * u}}}};
//   }
//   bool admissible(state_type const& s) const override {
//     auto const r = s[0];
//     auto const e = s[2];
//     return (r > 0 && e > 0);
//   }
//   state_type wave_speeds(state_type const& s) const override {
//     auto const p = pressure(s);
//     auto const r = s[0];
//     auto const u = s[1] / r;
//     auto const c = std::sqrt(m_params.gamma * p / r);
//     return state_type{u - c, u, u + c};
//   }
//   state_type prim_to_cons(state_type const& s) const override {
//     auto const& [r, u, p] = s;
//     return state_type{r, r * u, this->energy(s)};
//   }
//   state_type cons_to_prim(state_type const& s) const override {
//     auto const& [r, j, e] = s;
//     return state_type{r, j / r, this->pressure(s)};
//   }

// protected:
//   virtual value_type energy(state_type const& prim) const = 0;

//   Params m_params;
// };

template<std::size_t NPS = 0>
class EulerBase : public SystemBase<fivo::state<double, 3 + NPS>>,
                  public HasDensity<fivo::state<double, 3 + NPS>>,
                  public HasVelocity<fivo::state<double, 3 + NPS>>,
                  public HasPressure<fivo::state<double, 3 + NPS>> {
public:
  using base_type = SystemBase<fivo::state<double, 3 + NPS>>;
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

  static constexpr std::size_t num_passive_scalars() { return NPS; }

  struct Params { value_type gamma; };

  EulerBase(Mesh const& mesh,
            std::shared_ptr<BC> const& left_bc,
            std::shared_ptr<BC> const& right_bc,
            Params const& params)
    : base_type(mesh, left_bc, right_bc), m_params(params) {}

  value_type density(state_type const& s) const override { return s[0]; }
  value_type velocity(state_type const& s) const override { return s[1] / s[0]; }
  flux_type flux(state_type const& s) const override {
    auto const p = this->pressure(s);
    auto const& r = s[0];
    auto const& j = s[1];
    auto const& e = s[2];
    auto const u = j / r;
    flux_type fl;
    fl(0, 0) = j;
    fl(1, 0) = r * u * u + p;
    fl(2, 0) = (e + p) * u;
    for (std::size_t i = 0; i < NPS; ++i) { fl(3 + i, 0) = u * s[3 + i]; }
    return fl;
  }
  bool admissible(state_type const& s) const override {
    auto const r = s[0];
    auto const e = s[2];
    return (r > 0 && e > 0);
  }
  state_type wave_speeds(state_type const& s) const override {
    auto const p = this->pressure(s);
    auto const r = s[0];
    auto const u = s[1] / r;
    auto const c = std::sqrt(m_params.gamma * p / r);
    state_type ws = {u - c, u, u + c};
    for (std::size_t i = 0; i < NPS; ++i) { ws[3 + i] = u; }
    return ws;
  }
  state_type prim_to_cons(state_type const& s) const override {
    auto const& r = s[0];
    auto const& u = s[1];
    state_type cons = {r, r * u, energy(s)};
    for (std::size_t i = 0; i < NPS; ++i) { cons[3 + i] = r * s[3 + i]; }
    return cons;
  }
  state_type cons_to_prim(state_type const& s) const override {
    auto const& r = s[0];
    auto const& j = s[1];
    state_type prim = {r, j / r, this->pressure(s)};
    for (std::size_t i = 0; i < NPS; ++i) { prim[3 + i] = s[3 + i] / r; }
    return prim;
  }

protected:
  virtual value_type energy(state_type const& prim) const = 0;

  Params m_params;
};

} // namespace system

} // namespace fivo

#endif // FIVO_EULER_BASE_H
