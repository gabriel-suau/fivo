#ifndef FIVO_MATRIX_H
#define FIVO_MATRIX_H

#include <fivo/array.h>

#include <type_traits>
#include <tuple>

namespace fivo {

/** MATRIX HELPER */
template<typename T, std::size_t M, std::size_t N>
class static_matrix {
  static_assert(std::is_arithmetic<T>::value,
                "fivo::static_matrix : element type must be an arithmetic type.");
public:
  // Typedefs
  using value_type      = T;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference       = value_type&;
  using const_reference = value_type const&;
  using pointer         = value_type*;
  using const_pointer   = value_type const*;
  using iterator        = value_type*;
  using const_iterator  = value_type const*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
  T m_data[M * N];
  inline constexpr size_type index(size_type i, size_type j) const noexcept { return i * N + j; }

public:
  // Ctors/Dtor/Assignment
  static_matrix() = default;
  static_matrix(static_matrix const&) = default;
  static_matrix(static_matrix&&) noexcept = default;
  static_matrix& operator=(static_matrix const&) = default;
  static_matrix& operator=(static_matrix&&) noexcept = default;
  ~static_matrix() = default;

  explicit constexpr static_matrix(array<array<T, N>, M> const& init) {
    for (size_type i = 0; i < M; ++i)
      for (size_type j = 0; j < N; ++j)
        m_data[index(i, j)] = init[i][j];
  }

  // Element access
  constexpr const_reference operator()(size_type i, size_type j) const noexcept {
    return m_data[index(i, j)];
  }
  constexpr reference operator()(size_type i, size_type j) noexcept {
    return m_data[index(i, j)];
  }
  constexpr const_reference at(size_type i, size_type j) const {
    if (i >= rows() || j >= cols())
      throw std::out_of_range("fivo::static_matrix: Index out of range.");
    return m_data[index(i, j)];
  }
  constexpr reference at(size_type i, size_type j) {
    if (i >= rows() || j >= cols())
      throw std::out_of_range("fivo::static_matrix: Index out of range.");
    return m_data[index(i, j)];
  }
  constexpr const_reference front() const noexcept { return m_data[0]; }
  constexpr       reference front()       noexcept { return m_data[0]; }
  constexpr const_reference back()  const noexcept { return m_data[M * N - 1]; }
  constexpr       reference back()        noexcept { return m_data[M * N - 1]; }
  constexpr const_pointer   data()  const noexcept { return &m_data[0]; }
  constexpr       pointer   data()        noexcept { return &m_data[0]; }

  // Iterators
  constexpr       iterator begin()        noexcept { return       iterator(data()); }
  constexpr const_iterator begin()  const noexcept { return const_iterator(data()); }
  constexpr const_iterator cbegin() const noexcept { return const_iterator(data()); }
  constexpr       iterator end()          noexcept { return       iterator(data() + M * N); }
  constexpr const_iterator end()    const noexcept { return const_iterator(data() + M * N); }
  constexpr const_iterator cend()   const noexcept { return const_iterator(data() + M * N); }
  constexpr       reverse_iterator rbegin()        noexcept { return       reverse_iterator(end()); }
  constexpr const_reverse_iterator rbegin()  const noexcept { return const_reverse_iterator(end()); }
  constexpr const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(end()); }
  constexpr       reverse_iterator rend()          noexcept { return       reverse_iterator(begin() + M * N); }
  constexpr const_reverse_iterator rend()    const noexcept { return const_reverse_iterator(begin() + M * N); }
  constexpr const_reverse_iterator crend()   const noexcept { return const_reverse_iterator(begin() + M * N); }

  // Capacity
  constexpr bool empty() const noexcept { return size() == 0; }
  constexpr size_type size() const noexcept { return M * N; }
  constexpr size_type max_size() const noexcept { return size(); }
  constexpr size_type rows() const noexcept { return M; }
  constexpr size_type cols() const noexcept { return N; }

  // Other
  constexpr void fill(T const& value) {
    for (size_type i = 0; i < M * N; ++i) { m_data[i] = value; }
  }
  constexpr void swap(static_matrix& other) noexcept(std::is_nothrow_swappable_v<T>) {
    static_matrix tmp = other;
    other = *this;
    *this = tmp;
  }

  // Arithmetic operators
  constexpr inline auto& operator+=(static_matrix<T, M, N> const& rhs) {
    for (size_type i = 0; i < M * N; ++i) m_data[i] += rhs.m_data[i];
    return *this;
  }
  constexpr inline auto& operator-=(static_matrix<T, M, N> const& rhs) {
    for (size_type i = 0; i < M * N; ++i) m_data[i] -= rhs.m_data[i];
    return *this;
  }
  template<typename Scalar,
           std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
  constexpr inline auto& operator*=(Scalar const& factor) {
    for (size_type i = 0; i < M * N; ++i) m_data[i] *= factor;
    return *this;
  }
  template<typename Scalar,
           std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
  constexpr inline auto& operator/=(Scalar const& factor) {
    for (size_type i = 0; i < M * N; ++i) m_data[i] /= factor;
    return *this;
  }
};

template<typename T, std::size_t M, std::size_t N>
constexpr inline static_matrix<T, M, N> operator+(static_matrix<T, M, N> const& in) {
  return in;
}

template<typename T, std::size_t M, std::size_t N>
constexpr inline static_matrix<T, M, N> operator-(static_matrix<T, M, N> const& in) {
  static_matrix<T, M, N> out;
  for (typename static_matrix<T, M, N>::size_type i = 0; i < N; ++i) out[i] = -in[i];
  return out;
}

template<typename T, std::size_t M, std::size_t N>
constexpr inline static_matrix<T, M, N> operator+(static_matrix<T, M, N> lhs, static_matrix<T, M, N> const& rhs) {
  lhs += rhs;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N>
constexpr inline static_matrix<T, M, N> operator-(static_matrix<T, M, N> lhs, static_matrix<T, M, N> const& rhs) {
  lhs -= rhs;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N, typename Scalar,
         std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
constexpr inline static_matrix<T, M, N> operator*(static_matrix<T, M, N> lhs, Scalar const& fac) {
  lhs *= fac;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N, typename Scalar,
         std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
constexpr inline static_matrix<T, M, N> operator*(Scalar const& fac, static_matrix<T, M, N> lhs) {
  lhs *= fac;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N, typename Scalar>
constexpr inline static_matrix<T, M, N> operator/(static_matrix<T, M, N> lhs, Scalar const& fac) {
  lhs /= fac;
  return lhs;
}

// Matvec product
template<typename T, std::size_t M, std::size_t N>
constexpr inline array<T, M> operator*(static_matrix<T, M, N> const& lhs, array<T, N> const& rhs) {
  array<T, M> res;
  for (std::size_t i = 0; i < M; ++i) {
    res[i] = 0;
    for (std::size_t j = 0; j < N; ++j) {
      res[i] += lhs(i, j) * rhs[j];
    }
  }
  return res;
}

// Equality operators
template<typename T, std::size_t M, std::size_t N>
constexpr inline bool operator==(static_matrix<T, M, N> const& lhs,
                                 static_matrix<T, M, N> const& rhs) {
  for (std::size_t i = 0; i < M; ++i)
    for (std::size_t j = 0; j < N; ++j)
      if (lhs(i, j) != rhs(i, j)) return false;
  return true;
}
template<typename T, std::size_t M, std::size_t N>
constexpr inline bool operator!=(static_matrix<T, M, N> const& lhs,
                                 static_matrix<T, M, N> const& rhs) {
  return !(lhs == rhs);
}
template<typename T, std::size_t M, std::size_t N>
constexpr inline bool operator<(static_matrix<T, M, N> const& lhs,
                                static_matrix<T, M, N> const& rhs) {
  return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}
template<typename T, std::size_t M, std::size_t N>
constexpr inline bool operator>(static_matrix<T, M, N> const& lhs,
                                static_matrix<T, M, N> const& rhs) {
  return rhs < lhs;
}
template<typename T, std::size_t M, std::size_t N>
constexpr inline bool operator<=(static_matrix<T, M, N> const& lhs,
                                static_matrix<T, M, N> const& rhs) {
  return !(lhs > rhs);
}
template<typename T, std::size_t M, std::size_t N>
constexpr inline bool operator>=(static_matrix<T, M, N> const& lhs,
                                static_matrix<T, M, N> const& rhs) {
  return !(lhs < rhs);
}

} // namespace fivo

namespace std {

template<typename T, std::size_t M, std::size_t N>
void swap(::fivo::static_matrix<T, M, N>& lhs, ::fivo::static_matrix<T, M, N>& rhs)
  noexcept(noexcept(lhs.swap(rhs))) {
  lhs.swap(rhs);
}

} // namespace std

#endif // FIVO_MATRIX_H
