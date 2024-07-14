#ifndef FIVO_NUMFLUX_H
#define FIVO_NUMFLUX_H

#include <algorithm>
#include "system.h"
#include "mesh.h"
#include "traits.h"

namespace fivo {

namespace flux {

/**
 * \brief Base class template for numerical fluxes.
 *
 * \tparam Derived Type of the derived flux class (for CRTP)
 */
template<typename Derived>
struct NumericalFlux {

  /**
   * \brief Compute the numerical flux given a system and two states.
   *
   * \tparam System Type of the system
   *
   * \param[in] sys    System object for which to compute the numerical flux
   * \param[in] left   Left cell state
   * \param[in] right  Right cell state
   * \return Numerical flux computed using the \tparam Derived implementation
   *
   * The numerical flux is computed using the input system \a sys and the two
   * input states \a left and \a right.
   */
  template<typename System>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const {
    return static_cast<Derived const*>(this)->compute(sys, left, right);
  }

  /**
   * \brief Compute the numerical fluxes given a system, a global state and two
   * boundary conditions.
   *
   * \tparam System Type of the system
   *
   * \param[in] sys      System object for which to compute the numerical flux
   * \param[in] left_bc  State of the left boundary ghost-cell
   * \param[in] right_bc State of the right boundary ghost-cell
   * \param[in] state    Global state inside the domain
   * \return Numerical fluxes on the edges of the mesh contained in \a sys
   * computed using the \tparam Derived implementation.
   *
   * The numerical flux is also computed on the boundary edges using the
   * ghost-cell values \a left_bc and \a right_bc.
   */
  template<typename System>
  auto gcompute(System const& sys,
               typename System::state_type const& left_bc,
               typename System::state_type const& right_bc,
               typename System::global_state_type const& state) const {
    typename System::global_state_type numflux(state.size() + 1);
    numflux[0] = compute(sys, left_bc, state.front());
    for (std::size_t i = 1; i < numflux.size() - 1; ++i)
      numflux[i] = compute(sys, state[i-1], state[i]);
    numflux.back() = compute(sys, state.back(), right_bc);
    return numflux;
  }

  /**
   * \brief Add the numerical flux contribution
   *
   * \tparam System Type of the system
   *
   * \param[in] sys      System object (used only to ensure type correctness)
   * \param[in] numflux  Numerical fluxes on the edges (computed with tcompute)
   * \param[in] mesh     Mesh object
   * \param[out] source  Global source inside the domain
   */
  template<typename System>
  auto contribute(System const& sys, Mesh const& mesh,
                  typename System::global_state_type const& numflux,
                  typename System::global_state_type& source) const {
    source[0] += numflux[0];
    for (std::size_t i = 1; i < source.size(); ++i) {
      source[i - 1] -= numflux[i];
      source[i] += numflux[i];
    }
    source.back() -= numflux.back();
    source /= mesh.dx();
  }
};

/**
 * \brief Upwind numerical flux (only for \ref LinearAdvection).
 *
 * The results should match exactly those computed with \ref Godunov
 */
struct Upwind : NumericalFlux<Upwind> {
  static constexpr char const* name() { return "upwind"; }
  auto compute(system::LinearAdvection const& sys,
               typename system::LinearAdvection::state_type const& left,
               typename system::LinearAdvection::state_type const& right) const {
    auto const v = sys.velocity();
    return (v > 0) ? sys.flux(left) : sys.flux(right);
  }
};

/**
 * \brief Godunov numerical flux.
 */
struct Godunov : NumericalFlux<Godunov> {
  static constexpr char const* name() { return "godunov"; }
  template<typename System>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const {
    static_assert(traits::has_riemann_solver<System>::value,
                  "fivo::flux::Godunov::compute : input system goes not have an exact Riemann solver.");
    auto exact = sys.solve_riemann(left, right);
    return sys.flux(exact(0));
  }
};

/**
 * \brief Rusanov numerical flux.
 */
struct Rusanov : NumericalFlux<Rusanov> {
  static constexpr char const* name() { return "rusanov"; }
  template<typename System>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const {
    auto const fl = sys.flux(left);
    auto const fr = sys.flux(right);
    auto const wsl = sys.wave_speeds(left);
    auto const wsr = sys.wave_speeds(right);
    auto comp = [] (auto const& l, auto const& r) { return std::abs(l) < std::abs(r); };
    auto const lmax = *std::max_element(wsl.begin(), wsl.end(), comp);
    auto const rmax = *std::max_element(wsr.begin(), wsr.end(), comp);
    auto const c = std::max(std::abs(lmax), std::abs(rmax));
    return 0.5 * (fr + fl - c * (right - left));
  }
};

/**
 * \brief HLL numerical flux.
 */
struct HLL : NumericalFlux<HLL> {
  static constexpr char const* name() { return "hll"; }
  template<typename System>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const {
    auto const fl = sys.flux(left);
    auto const fr = sys.flux(right);
    auto const wsl = sys.wave_speeds(left);
    auto const wsr = sys.wave_speeds(right);
    auto const lmin = *std::min_element(wsl.begin(), wsl.end());
    auto const rmin = *std::min_element(wsr.begin(), wsr.end());
    auto const c1 = std::min(lmin, rmin);
    if (0 <= c1) return fl;
    auto const lmax = *std::max_element(wsl.begin(), wsl.end());
    auto const rmax = *std::max_element(wsr.begin(), wsr.end());
    auto const c2 = std::max(lmax, rmax);
    if (c2 <= 0) return fr;
    return (c2 * fl - c1 * fr + c1 * c2 * (right - left)) / (c2 - c1);
  }
};

/**
 * \brief HLLC numerical flux.
 */
struct HLLC : NumericalFlux<HLLC> {
  static constexpr char const* name() { return "hllc"; }
  template<typename System,
           std::enable_if_t<traits::is_derived<System, system::Euler>::value, bool> = false>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const {
    auto const fl = sys.flux(left);
    auto const fr = sys.flux(right);
    auto const wsl = sys.wave_speeds(left);
    auto const wsr = sys.wave_speeds(right);
    auto const lmin = *std::min_element(wsl.begin(), wsl.end());
    auto const rmin = *std::min_element(wsr.begin(), wsr.end());
    auto const cl = std::min(lmin, rmin);
    if (0 <= cl) return fl;
    auto const lmax = *std::max_element(wsl.begin(), wsl.end());
    auto const rmax = *std::max_element(wsr.begin(), wsr.end());
    auto const cr = std::max(lmax, rmax);
    if (cr <= 0) return fr;

    // Compute intermediate state
    auto const& [rl, ul, pl] = sys.cons_to_prim(left);
    auto const& [rr, ur, pr] = sys.cons_to_prim(right);
    auto const cstar = (pr - pl + rl * ul * (cl - ul) - rr * ur * (cr - ur))
      / (rl * (cl - ul) - rr * (cr - ur));
    auto const dstar = typename System::state_type{0, 1, cstar};
    auto const plstar = pl + rl * (cl - ul) * (cstar - ul);
    auto const lstar  = (cl * left  - fl + plstar * dstar) / (cl - cstar);
    auto const flstar = fl + cl * (lstar - left);
    if (0 <= cstar) return flstar;
    auto const prstar = pr + rr * (cr - ur) * (cstar - ur);
    auto const rstar  = (cr * right - fr + prstar * dstar) / (cr - cstar);
    auto const frstar = fr + cr * (rstar - right);
    return frstar;
  }

  template<template<std::size_t> typename System, std::size_t NPS,
           std::enable_if_t<traits::is_derived<System<NPS>, system::EulerP<NPS>>::value,
                            bool> = false>
  auto compute(System<NPS> const& sys,
               typename System<NPS>::state_type const& left,
               typename System<NPS>::state_type const& right) const {
    auto const fl = sys.flux(left);
    auto const fr = sys.flux(right);
    auto const wsl = sys.wave_speeds(left);
    auto const wsr = sys.wave_speeds(right);
    auto const lmin = *std::min_element(wsl.begin(), wsl.end());
    auto const rmin = *std::min_element(wsr.begin(), wsr.end());
    auto const cl = std::min(lmin, rmin);
    if (0 <= cl) return fl;
    auto const lmax = *std::max_element(wsl.begin(), wsl.end());
    auto const rmax = *std::max_element(wsr.begin(), wsr.end());
    auto const cr = std::max(lmax, rmax);
    if (cr <= 0) return fr;

    // Compute intermediate state
    auto const lprim = sys.cons_to_prim(left);
    auto const rprim = sys.cons_to_prim(right);
    auto const& rl = lprim[0];
    auto const& ul = lprim[1];
    auto const& pl = lprim[2];
    auto const& rr = rprim[0];
    auto const& ur = rprim[1];
    auto const& pr = rprim[2];
    auto const cstar = (pr - pl + rl * ul * (cl - ul) - rr * ur * (cr - ur))
      / (rl * (cl - ul) - rr * (cr - ur));
    auto dstar = typename System<NPS>::state_type{0, 1, cstar};
    for (std::size_t i = 0; i < NPS; ++i) { dstar[3 + i] = 0.; }
    auto const plstar = pl + rl * (cl - ul) * (cstar - ul);
    auto const lstar  = (cl * left  - fl + plstar * dstar) / (cl - cstar);
    auto const flstar = fl + cl * (lstar - left);
    if (0 <= cstar) return flstar;
    auto const prstar = pr + rr * (cr - ur) * (cstar - ur);
    auto const rstar  = (cr * right - fr + prstar * dstar) / (cr - cstar);
    auto const frstar = fr + cr * (rstar - right);
    return frstar;
  }
};
};

} // namespace flux

} // namespace fivo

#endif // FIVO_NUMFLUX_H
