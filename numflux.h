#ifndef FIVO_NUMFLUX_H
#define FIVO_NUMFLUX_H

#include <algorithm>
#include "system.h"
#include "mesh.h"

namespace fivo {

template<typename Derived>
struct NumericalFlux {
  template<typename System>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const {
    return static_cast<Derived const*>(this)->compute(sys, left, right);
  }
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
};

struct Rusanov : NumericalFlux<Rusanov> {
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

struct HLL : NumericalFlux<HLL> {
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

struct HLLC : NumericalFlux<HLLC> {
  template<typename System>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const;

  auto compute(IdealGasEuler const& sys,
               typename IdealGasEuler::state_type const& left,
               typename IdealGasEuler::state_type const& right) const {
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
    auto const rl = sys.density(left);
    auto const rr = sys.density(right);
    auto const pl = sys.pressure(left);
    auto const pr = sys.pressure(right);
    auto const ul = sys.velocity(left);
    auto const ur = sys.velocity(right);
    auto const cstar = (pr - pl + rl * ul * (cl - ul) - rr * ur * (cr - ur))
      / (rl * (cl - ul) - rr * (cr - ur));
    auto const dstar = typename IdealGasEuler::state_type{0, 1, cstar};
    auto const plstar = pl + rl * (cl - ul) * (cstar - ul);
    auto const ulstar = (cl * left  - fl + plstar * dstar) / (cl - cstar);
    auto const flstar = fl + cl * (ulstar - left);
    if (0 <= cstar) return flstar;
    auto const prstar = pr + rr * (cr - ur) * (cstar - ur);
    auto const urstar = (cr * right - fr + prstar * dstar) / (cr - cstar);
    auto const frstar = fr + cr * (urstar - right);
    return frstar;
  }
};

} // namespace fivo

#endif // FIVO_NUMFLUX_H
