#ifndef FIVO_BASIS_H
#define FIVO_BASIS_H

#include "matrix.h"
#include <vector>
#include <algorithm>
#include <cmath>

namespace fivo {

template<typename Scalar>
struct Polynomial {
  static_assert(std::is_arithmetic<Scalar>::value,
                "fivo::Polynomial : Scalar type must be an arithmetic type.");
  using type = Polynomial<Scalar>;
  using scalar_type = Scalar;

  Polynomial() : m_coef(1, 0.) {}
  Polynomial(scalar_type const& c) : m_coef{c} {}
  explicit Polynomial(std::size_t degree) : m_coef(degree + 1, 0.) {}
  explicit Polynomial(std::vector<Scalar> const& coef) : m_coef(coef) {}
  explicit Polynomial(std::initializer_list<Scalar> coef) : m_coef(coef) {}

  auto degree() const { return m_coef.size() - 1; }
  auto const& coefficients() const { return m_coef; }

  type derivative() const {
    if (degree() == 0) return type();
    auto deriv = type(degree() - 1);
    for (std::size_t i = 1; i <= degree(); ++i) { deriv.m_coef[i - 1] = i * m_coef[i]; }
    return deriv;
  }

  type primitive() const {
    auto prim = type(degree() + 1);
    prim.m_coef[0] = 0.;
    for (std::size_t i = 1; i <= prim.degree(); ++i) { prim.m_coef[i] = m_coef[i - 1] / i; }
    return prim;
  }

  auto operator()(scalar_type const& p) const {
    scalar_type eval = m_coef[degree()];
    for (int i = degree() - 1; i >= 0; --i) { eval = p * eval + m_coef[i]; }
    return eval;
  }

  /* ARITHMETIC OPERATIONS */
  auto operator+() const { return *this; }
  auto operator-() const {
    auto copy = type(*this);
    std::transform(m_coef.begin(), m_coef.end(), copy.m_coef.begin(),
                   [] (auto const& c) { return -c; });
    return copy;
  }
  auto& operator+=(type const& rhs) {
    auto const lim = rhs.degree() + 1;
    if (degree() < rhs.degree()) m_coef.resize(lim, 0.);
    for (int i = 0; i < lim; ++i) m_coef[i] += rhs.m_coef[i];
    return *this;
  }
  auto& operator-=(type const& rhs) {
    auto const lim = rhs.degree() + 1;
    if (degree() < rhs.degree()) m_coef.resize(lim, 0.);
    for (int i = 0; i < lim; ++i) m_coef[i] -= rhs.m_coef[i];
    return *this;
  }
  auto& operator*=(type const& rhs) {
    auto old = type(*this);
    auto const newdeg = degree() + rhs.degree();
    m_coef.resize(newdeg + 1);
    std::fill(m_coef.begin(), m_coef.end(), 0.);
    for (int i = 0; i <= old.degree(); ++i)
      for (int j = 0; j <= rhs.degree(); ++j)
        m_coef[i + j] += old.m_coef[i] * rhs.m_coef[j];
    return *this;
  }
  template<typename T,
           std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
  auto& operator*=(T const& s) {
    std::for_each(m_coef.begin(), m_coef.end(), [&] (auto& c) { c *= s; });
    return *this;
  }
  template<typename T,
           std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
  auto& operator/=(T const& s) {
    std::for_each(m_coef.begin(), m_coef.end(), [&] (auto& c) { c /= s; });
    return *this;
  }
private:
  std::vector<Scalar> m_coef;
};

template<typename Scalar>
auto operator+(Polynomial<Scalar> lhs, Polynomial<Scalar> const& rhs) {
  lhs += rhs;
  return lhs;
}
template<typename Scalar>
auto operator-(Polynomial<Scalar> lhs, Polynomial<Scalar> const& rhs) {
  lhs -= rhs;
  return lhs;
}
template<typename Scalar>
auto operator*(Polynomial<Scalar> lhs, Polynomial<Scalar> const& rhs) {
  lhs *= rhs;
  return lhs;
}
template<typename Scalar, typename T>
auto operator*(Polynomial<Scalar> lhs, T const& s) {
  lhs *= s;
  return lhs;
}
template<typename Scalar, typename T>
auto operator*(T const& s, Polynomial<Scalar> lhs) {
  lhs *= s;
  return lhs;
}
template<typename Scalar, typename T>
auto operator/(Polynomial<Scalar> lhs, T const& s) {
  lhs /= s;
  return lhs;
}

template<typename Poly>
struct Integrator {
  using scalar_type = typename Poly::scalar_type;
  static_assert(std::is_arithmetic<scalar_type>::value,
                "fivo::Integrator : Poly::scalar_type must be an arithmetic type.");
  scalar_type integrate(Poly const& poly, std::pair<scalar_type, scalar_type> const& range) {
    auto const prim = poly.primitive();
    return prim(range.second) - prim(range.first);
  }
};

template<typename Poly, typename T>
typename Poly::scalar_type integrate(Poly const& poly, T const& a, T const& b) {
  static_assert(std::is_arithmetic<T>::value,
                "fivo:integrate : T must be an arithmetic type.");
  auto const prim = poly.primitive();
  return prim(b) - prim(a);
}

template<typename Poly, typename T>
typename Poly::scalar_type integrate(Poly const& poly, std::pair<T, T> const& range) {
  static_assert(std::is_arithmetic<T>::value,
                "fivo:integrate : range must be a pair of arithmetic type.");
  auto const prim = poly.primitive();
  return prim(range.second) - prim(range.first);
}

struct Basis {
  using func_type = Polynomial<double>;
  explicit Basis(int order) : m_func(order + 1) {}
  auto size() const { return m_func.size(); }
  auto order() const { return size() - 1; }
  auto const& functions() const { return m_func; }
protected:
  std::vector<func_type> m_func;
  matrix<double> mass, stiffness;
};

struct LegendreBasis : Basis {
  explicit LegendreBasis(int order, bool normalise = true) : Basis(order) {
    build();
    if (normalise) this->normalise();
  }
private:
  void build() {
    m_func[0] = 1;
    if (order() == 0) return;
    m_func[1] = func_type({0, 1});
    for (std::size_t i = 2; i <= order(); ++i) {
      m_func[i] = (m_func[i - 1] * (2 * i - 1) * func_type({0, 1})
                   - m_func[i - 2] * (i - 1)) / i;
    }
  }
  void normalise() {
    auto dot = [&] (auto const& lhs, auto const& rhs) { return integrate(lhs * rhs, -1, 1); };
    for (auto& f : m_func) { f /= std::sqrt(dot(f, f)); }
  }
};

// struct LegendreBasis : Basis {
  
// };

} // namespace fivo

#endif // FIVO_BASIS_H
