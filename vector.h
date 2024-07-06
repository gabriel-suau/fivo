#ifndef FIVO_VECTOR_H
#define FIVO_VECTOR_H

#include <vector>

namespace fivo {

template<typename T>
using vector = std::vector<T>;

template<typename State>
auto& operator+=(vector<State>& lhs, vector<State> const& rhs) {
  for (typename vector<State>::size_type i = 0; i < lhs.size(); ++i) lhs[i] += rhs[i];
  return lhs;
}

template<typename State>
auto& operator-=(vector<State>& lhs, vector<State> const& rhs) {
  for (typename vector<State>::size_type i = 0; i < lhs.size(); ++i) lhs[i] -= rhs[i];
  return lhs;
}

template<typename State, typename Scalar>
auto& operator*=(vector<State>& lhs, Scalar const& factor) {
  for (typename vector<State>::size_type i = 0; i < lhs.size(); ++i) lhs[i] *= factor;
  return lhs;
}

template<typename State, typename Scalar>
auto& operator/=(vector<State>& lhs, Scalar const& factor) {
  for (typename vector<State>::size_type i = 0; i < lhs.size(); ++i) lhs[i] /= factor;
  return lhs;
}

template<typename State>
auto operator+(vector<State> const& in) {
  return in;
}

template<typename State>
auto operator-(vector<State> const& in) {
  vector<State> out(in.size());
  for (typename vector<State>::size_type i = 0; i < in.size(); ++i) out[i] = -in[i];
  return out;
}

template<typename State>
auto operator+(vector<State> lhs, vector<State> const& rhs) {
  lhs += rhs;
  return lhs;
}

template<typename State>
auto operator-(vector<State> lhs, vector<State> const& rhs) {
  lhs -= rhs;
  return lhs;
}

template<typename State, typename Scalar>
auto operator*(vector<State> lhs, Scalar const& fac) {
  lhs *= fac;
  return lhs;
}

template<typename State, typename Scalar>
auto operator*(Scalar const& fac, vector<State> rhs) {
  return rhs * fac;
}

template<typename State, typename Scalar>
auto operator/(vector<State> lhs, Scalar const& fac) {
  lhs /= fac;
  return lhs;
}

} // namespace fivo

#endif // FIVO_VECTOR_H
