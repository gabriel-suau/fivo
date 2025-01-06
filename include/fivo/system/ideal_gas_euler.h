#ifndef FIVO_IDEAL_GAS_EULER_H
#define FIVO_IDEAL_GAS_EULER_H

#include <fivo/system/euler_base.h>
#include <iostream>

namespace fivo {

namespace system {

// class IdealGasEuler : public EulerBase,
//                       public HasRiemannSolver<IdealGasEuler,
//                                               typename EulerBase::state_type> {
// public:
//   using base_type = EulerBase;
//   using state_type = typename base_type::state_type;
//   using flux_type = typename base_type::flux_type;
//   using global_state_type = typename base_type::global_state_type;
//   using value_type = typename base_type::value_type;

//   using Params = typename base_type::Params;
//   using BC = typename EulerBase::BC;

//   using EulerBase::EulerBase;
//   // IdealGasEuler(Mesh const& mesh,
//   //               std::shared_ptr<BC> const& left_bc,
//   //               std::shared_ptr<BC> const& right_bc,
//   //               Params const& params)
//   //   : base_type(mesh, left_bc, right_bc, {}), m_params(params) {}

//   // auto const& get_params() const { return m_params; }
//   // void set_params(Params const& params) { m_params = params; }

//   value_type pressure(state_type const& s) const override {
//     auto const& [r, j, e] = s;
//     return (this->m_params.gamma - 1) * (e - 0.5 * j * j / r);
//   }

//   auto solve_riemann(state_type const& left, state_type const& right) const {
//     auto const& gamma = this->m_params.gamma;
//     auto const& rl = left[0];
//     auto const& rr = right[0];
//     auto const pl = pressure(left);
//     auto const pr = pressure(right);
//     auto const ul = velocity(left);
//     auto const ur = velocity(right);
//     auto const cl = std::sqrt(gamma * pl / rl);
//     auto const cr = std::sqrt(gamma * pr / rr);

//     // Find pstar and ustar in the star region
//     auto const Al = 2. / (rl * (gamma + 1));
//     auto const Ar = 2. / (rr * (gamma + 1));
//     auto const Bl = pl * (gamma - 1) / (gamma + 1);
//     auto const Br = pr * (gamma - 1) / (gamma + 1);
//     auto const du = ur - ul;

//     auto const fk =
//       [=] (auto const p, auto const pk, auto const ck, auto const Ak, auto const Bk) {
//         if (p > pk) return (p - pk) * std::sqrt(Ak / (p + Bk));
//         else return 2 * ck / (gamma - 1) * (std::pow(p / pk, (gamma - 1) / (2 * gamma)) - 1);
//       };
//     auto const dfk =
//       [=] (auto const p, auto const pk, auto const rk,
//           auto const ck, auto const Ak, auto const Bk) {
//         if (p > pk) return std::sqrt(Ak / (Bk + p)) * (1 - (p - pk) / (2 * (Bk + p)));
//         else return std::pow(p / pk, - (gamma + 1) / (2 * gamma)) / (rk * ck);
//       };
//     auto const fl = [&] (auto const p) { return fk(p, pl, cl, Al, Bl); };
//     auto const fr = [&] (auto const p) { return fk(p, pr, cr, Ar, Br); };
//     auto const dfl = [&] (auto const p) { return dfk(p, pl, rl, cl, Al, Bl); };
//     auto const dfr = [&] (auto const p) { return dfk(p, pr, rr, cr, Ar, Br); };
//     auto const f = [&] (auto const p) { return fl(p) + fr(p) + du; };
//     auto const df = [&] (auto const p) { return dfl(p) + dfr(p); };

//     auto const tol = 1e-10;
//     auto const ppv = 0.5 * (pl + pr) - 0.125 * (ur - ul) * (rl + rr) * (cl + cr);
//     auto const p0 = std::max(tol, ppv);

//     auto const pstar = math::newton_raphson(p0, f, df, tol);
//     auto const ustar = 0.5 * ((ul + ur) + fr(pstar) - fl(pstar));

//     /* FIND rlstar AND rrstar (density in the star region) */
//     auto const fac = (gamma - 1) / (gamma + 1);

//     // Solution left of the contact
//     auto const rlstar = (pstar > pl)
//       ? rl * (pstar / pl + fac) / ((pstar / pl) * fac + 1) // left shock
//       : rl * std::pow(pstar / pl, 1 / gamma);            // left rarefaction

//     // Solution right of the contact
//     auto const rrstar = (pstar > pr)
//       ? rr * (pstar / pr + fac) / ((pstar / pr) * fac + 1) // right shock
//       : rr * std::pow(pstar / pr, 1 / gamma);            // right rarefaction

//     // Left shock speed
//     auto const sl = ul - cl * std::sqrt(((gamma + 1) * pstar / pl + (gamma - 1)) / (2 * gamma));
//     // Left head and tail rarefaction speeds
//     auto const clstar = cl * std::pow(pstar / pl, (gamma - 1) / (2 * gamma));
//     auto const shl = ul - cl;
//     auto const stl = ustar - clstar;

//     // Right shock speed
//     auto const sr = ur + cr * std::sqrt(((gamma + 1) * pstar / pr + (gamma - 1)) / (2 * gamma));
//     // Right head and tail rarefaction speeds
//     auto const crstar = cr * std::pow(pstar / pr, (gamma - 1) / (2 * gamma));
//     auto const shr = ur + cr;
//     auto const str = ustar + crstar;

//     // Solution in the left star region\fan
//     auto const lstar =
//       [=] (value_type const&) { return this->prim_to_cons(state_type{rlstar, ustar, pstar}); };

//     // Solution in the right star region\fan
//     auto const rstar =
//       [=] (value_type const&) { return this->prim_to_cons(state_type{rrstar, ustar, pstar}); };

//     // Solution in the fan region for left rarefaction
//     auto const lfan =
//       [=] (value_type const& xt) {
//         auto const tmp = (2 + (gamma - 1) * (ul - xt) / cl) / (gamma + 1);
//         auto const rlfan = rl * std::pow(tmp, 2 / (gamma - 1));
//         auto const ulfan = 2 / (gamma + 1) * (cl + (gamma - 1) * ul / 2 + xt);
//         auto const plfan = pl * std::pow(tmp, 2 * gamma / (gamma - 1));
//         return this->prim_to_cons(state_type{rlfan, ulfan, plfan});
//       };

//     // Solution in the fan region for right rarefaction
//     auto const rfan =
//       [=] (value_type const& xt) {
//         auto const tmp = (2 - (gamma - 1) * (ur - xt) / cr) / (gamma + 1);
//         auto const rrfan = rr * std::pow(tmp, 2 / (gamma - 1));
//         auto const urfan = 2 / (gamma + 1) * (-cr + (gamma - 1) * ur / 2 + xt);
//         auto const prfan = pr * std::pow(tmp, 2 * gamma / (gamma - 1));
//         return this->prim_to_cons(state_type{rrfan, urfan, prfan});
//       };

//     auto const sol =
//       [=] (value_type const& xt) {
//         if (xt < ustar) { // Left of contact
//           if (pstar > pl) { return (xt < sl) ? left : lstar(xt); }
//           else { return (xt < shl) ? left : (xt > stl) ? lstar(xt) : lfan(xt); }
//         }
//         else { // Right of contact
//           if (pstar > pr) { return (xt > sr) ? right : rstar(xt);}
//           else { return (xt > shr) ? right : (xt < str) ? rstar(xt) : rfan(xt); }
//         }
//       };
//     return sol;
//   }

// private:
//   value_type energy(state_type const& prim) const override {
//     auto const& [r, u, p] = prim;
//     return p / (this->m_params.gamma - 1) + 0.5 * r * u * u;
//   }

//   // Params m_params;
// };

template<std::size_t NPS = 0>
class IdealGasEuler : public EulerBase<NPS>,
                      public HasRiemannSolver<IdealGasEuler<NPS>,
                                              typename EulerBase<NPS>::state_type> {
public:
  using base_type = EulerBase<NPS>;
  using state_type = typename base_type::state_type;
  using flux_type = typename base_type::flux_type;
  using global_state_type = typename base_type::global_state_type;
  using value_type = typename base_type::value_type;

  using Params = typename base_type::Params;
  using BC = typename base_type::BC;

  using EulerBase<NPS>::EulerBase;
  // IdealGasEulerP(Mesh const& mesh,
  //                std::shared_ptr<BC> const& left_bc,
  //                std::shared_ptr<BC> const& right_bc,
  //                Params const& params)
  //   : base_type(mesh, left_bc, right_bc, params) {}

  // auto const& get_params() const { return m_params; }
  // void set_params(Params const& params) { m_params = params; }

  value_type pressure(state_type const& s) const override {
    auto const& r = s[0];
    auto const& j = s[1];
    auto const& e = s[2];
    return (this->m_params.gamma - 1) * (e - 0.5 * j * j / r);
  }

  auto solve_riemann(state_type const& left, state_type const& right) const {
    auto const& gamma = this->m_params.gamma;
    auto const pleft = this->cons_to_prim(left);
    auto const pright = this->cons_to_prim(right);
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
        for (std::size_t i = 0; i < NPS; ++i) { prim[3 + i] = pleft[3 + i]; }
        return this->prim_to_cons(prim);
      };

    // Solution in the right star region\fan
    auto const rstar =
      [=] (value_type const&) {
        state_type prim = {rrstar, ustar, pstar};
        for (std::size_t i = 0; i < NPS; ++i) { prim[3 + i] = pright[3 + i]; }
        return this->prim_to_cons(prim);
      };

    // Solution in the fan region for left rarefaction
    auto const lfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 + (gamma - 1) * (ul - xt) / cl) / (gamma + 1);
        auto const rlfan = rl * std::pow(tmp, 2 / (gamma - 1));
        auto const ulfan = 2 / (gamma + 1) * (cl + (gamma - 1) * ul / 2 + xt);
        auto const plfan = pl * std::pow(tmp, 2 * gamma / (gamma - 1));
        state_type prim = {rlfan, ulfan, plfan};
        for (std::size_t i = 0; i < NPS; ++i) { prim[3 + i] = pleft[3 + i]; }
        return this->prim_to_cons(prim);
      };

    // Solution in the fan region for right rarefaction
    auto const rfan =
      [=] (value_type const& xt) {
        auto const tmp = (2 - (gamma - 1) * (ur - xt) / cr) / (gamma + 1);
        auto const rrfan = rr * std::pow(tmp, 2 / (gamma - 1));
        auto const urfan = 2 / (gamma + 1) * (-cr + (gamma - 1) * ur / 2 + xt);
        auto const prfan = pr * std::pow(tmp, 2 * gamma / (gamma - 1));
        state_type prim = {rrfan, urfan, prfan};
        for (std::size_t i = 0; i < NPS; ++i) { prim[3 + i] = pright[3 + i]; }
        return this->prim_to_cons(prim);
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

private:
  value_type energy(state_type const& prim) const override {
    auto const& r = prim[0];
    auto const& u = prim[1];
    auto const& p = prim[2];
    return p / (this->m_params.gamma - 1) + 0.5 * r * u * u;
  }

  // Params m_params;
};

} // namespace system

} // namespace fivo

#endif // FIVO_IDEAL_GAS_EULER_H
