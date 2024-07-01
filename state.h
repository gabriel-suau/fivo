#ifndef FIVO_STATE_H
#define FIVO_STATE_H

#include "array.h"
#include <vector>

namespace fivo {

/* LOCAL STATE */
template<typename T, std::size_t DIM>
using state = array<T, DIM>;

/* GLOBAL STATE */
template<typename State>
using global_state = std::vector<State>;

template<typename State>
auto& operator+=(fivo::global_state<State>& lhs, fivo::global_state<State> const& rhs) {
  for (int i = 0; i < lhs.size(); ++i) lhs[i] += rhs[i];
  return lhs;
}

template<typename State>
auto& operator-=(fivo::global_state<State>& lhs, fivo::global_state<State> const& rhs) {
  for (int i = 0; i < lhs.size(); ++i) lhs[i] -= rhs[i];
  return lhs;
}

template<typename State, typename Scalar>
auto& operator*=(fivo::global_state<State>& lhs, Scalar const& factor) {
  for (int i = 0; i < lhs.size(); ++i) lhs[i] *= factor;
  return lhs;
}

template<typename State, typename Scalar>
auto& operator/=(fivo::global_state<State>& lhs, Scalar const& factor) {
  for (int i = 0; i < lhs.size(); ++i) lhs[i] /= factor;
  return lhs;
}

template<typename State>
auto operator-(fivo::global_state<State> const& in) {
  fivo::global_state<State> out(in.size());
  for (int i = 0; i < in.size(); ++i) out[i] = -in[i];
  return out;
}

template<typename State>
auto operator+(fivo::global_state<State> const& in) {
  return in;
}

template<typename State>
auto operator+(fivo::global_state<State> lhs, fivo::global_state<State> const& rhs) {
  lhs += rhs;
  return lhs;
}

template<typename State>
auto operator-(fivo::global_state<State> lhs, fivo::global_state<State> const& rhs) {
  lhs -= rhs;
  return lhs;
}

template<typename State, typename Scalar>
auto operator*(fivo::global_state<State> lhs, Scalar const& fac) {
  lhs *= fac;
  return lhs;
}

template<typename State, typename Scalar>
auto operator*(Scalar const& fac, fivo::global_state<State> rhs) {
  return rhs * fac;
}

template<typename State, typename Scalar>
auto operator/(fivo::global_state<State> lhs, Scalar const& fac) {
  lhs /= fac;
  return lhs;
}

} // namespace fivo

#endif // FIVO_STATE_H
