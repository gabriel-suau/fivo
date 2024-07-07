#ifndef FIVO_MATRIX_H
#define FIVO_MATRIX_H

#include <type_traits>
#include <vector>
#include <stdexcept>

namespace fivo {

template<typename T>
struct matrix {
  static_assert(std::is_arithmetic<T>::value,
                "fivo::matrix : element type must be an arithmetic type.");

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

  // Constructors

  // Data
  size_type m_rows, m_cols;
  std::vector<T> storage;

  // Iterators
  constexpr const_iterator begin()  const noexcept { return const_iterator(data()); }
  constexpr iterator       begin()  noexcept       { return iterator(data()); }
  constexpr const_iterator cbegin() const noexcept { return const_iterator(data()); }
  constexpr const_iterator end()  const noexcept { return const_iterator(data() + m_rows * m_cols); }
  constexpr iterator       end()  noexcept       { return iterator(data() + m_rows * m_cols); }
  constexpr const_iterator cend() const noexcept { return const_iterator(data() + m_rows * m_cols); }

  // Access
  constexpr const_reference operator()(size_type i, size_type j) const noexcept { return storage[i * cols() + j]; }
  constexpr reference       operator()(size_type i, size_type j)       noexcept { return storage[i * cols() + j]; }

  constexpr const_reference front() const noexcept { return storage[0]; }
  constexpr reference       front()       noexcept { return storage[0]; }

  constexpr const_reference back() const noexcept { return storage[size() - 1]; }
  constexpr reference       back()       noexcept { return storage[size() - 1]; }

  constexpr const_pointer data() const noexcept { return &storage[0]; }
  constexpr pointer       data()       noexcept { return &storage[0]; }

  constexpr const_reference at(size_type i, size_type j) const {
    if (i >= m_rows) throw std::out_of_range("fivo::matrix::at : row index is out of bounds.");
    if (j >= m_cols) throw std::out_of_range("fivo::matrix::at : col index is out of bounds.");
    return storage[i * cols() + j];
  }
  constexpr reference at(size_type i, size_type j) {
    if (i >= m_rows) throw std::out_of_range("fivo::matrix::at : row index is out of bounds.");
    if (j >= m_cols) throw std::out_of_range("fivo::matrix::at : col index is out of bounds.");
    return storage[i * cols() + j];
  }

  // Capacity
  constexpr bool empty() const { return begin() == end(); }
  constexpr size_type rows() const noexcept { return m_rows; }
  constexpr size_type cols() const noexcept { return m_cols; }
  constexpr size_type size() const noexcept { return m_rows * m_cols; }
  constexpr size_type max_size() const noexcept { return size(); }
  constexpr void resize(size_type rows, size_type cols) {
    m_rows = rows;
    m_cols = cols;
    storage.resize(rows * cols);
  }
  constexpr void clear() noexcept { m_rows = 0; m_cols = 0; storage.clear(); }
  // Other
  constexpr void fill(T const& value) {
    for (size_type i = 0; i < size(); ++i) storage[i] = value;
  }

  constexpr void swap(matrix<T>& other) {
    matrix<T> tmp = other;
    other = *this;
    *this = tmp;
  }

  // Arithmetic operators
  auto& operator+=(matrix<T> const& rhs) {
    for (size_type i = 0; i < size(); ++i) storage[i] += rhs.storage[i];
    return *this;
  }
  auto& operator-=(matrix<T> const& rhs) {
    for (size_type i = 0; i < size(); ++i) storage[i] -= rhs.storage[i];
    return *this;
  }
  template<typename Scalar>
  auto& operator*=(Scalar const& factor) {
    for (size_type i = 0; i < size(); ++i) storage[i] *= factor;
    return *this;
  }
  template<typename Scalar>
  auto& operator/=(Scalar const& factor) {
    for (size_type i = 0; i < size(); ++i) storage[i] /= factor;
    return *this;
  }
};

/** MATRIX HELPER */
template<typename T, std::size_t M, std::size_t N>
struct static_matrix {
  static_assert(std::is_arithmetic<T>::value,
                "fivo::static_matrix : element type must be an arithmetic type.");

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

  // Data
  T storage[M * N];

  // Iterators
  constexpr const_iterator begin()  const noexcept { return const_iterator(data()); }
  constexpr iterator       begin()  noexcept       { return iterator(data()); }
  constexpr const_iterator cbegin() const noexcept { return const_iterator(data()); }
  constexpr const_iterator end()  const noexcept { return const_iterator(data() + M * N); }
  constexpr iterator       end()  noexcept       { return iterator(data() + M * N); }
  constexpr const_iterator cend() const noexcept { return const_iterator(data() + M * N); }

  // Access
  constexpr const_reference operator()(size_type i, size_type j) const noexcept { return storage[i * N + j]; }
  constexpr reference       operator()(size_type i, size_type j)       noexcept { return storage[i * N + j]; }

  constexpr const_reference front() const noexcept { return storage[0]; }
  constexpr reference       front()       noexcept { return storage[0]; }

  constexpr const_reference back() const noexcept { return storage[M * N - 1]; }
  constexpr reference       back()       noexcept { return storage[M * N - 1]; }

  constexpr const_pointer data() const noexcept { return &storage[0]; }
  constexpr pointer       data()       noexcept { return &storage[0]; }

  // Capacity
  constexpr size_type rows() const noexcept { return M; }
  constexpr size_type cols() const noexcept { return N; }
  constexpr size_type size() const noexcept { return M * N; }
  constexpr size_type max_size() const noexcept { return M * N; }

  // Other
  constexpr void fill(T const& value) {
    for (size_type i = 0; i < M * N; ++i) storage[i] = value;
  }
  constexpr void swap(static_matrix& other) {
    static_matrix tmp = other;
    other = *this;
    *this = tmp;
  }

  // Arithmetic operators
  auto& operator+=(static_matrix<T, M, N> const& rhs) {
    for (size_type i = 0; i < M * N; ++i) storage[i] += rhs[i];
    return *this;
  }
  auto& operator-=(static_matrix<T, M, N> const& rhs) {
    for (size_type i = 0; i < M * N; ++i) storage[i] -= rhs[i];
    return *this;
  }
  template<typename Scalar>
  auto& operator*=(Scalar const& factor) {
    for (size_type i = 0; i < M * N; ++i) storage[i] *= factor;
    return *this;
  }
  template<typename Scalar>
  auto& operator/=(Scalar const& factor) {
    for (size_type i = 0; i < M * N; ++i) storage[i] /= factor;
    return *this;
  }
};

template<typename T, std::size_t M, std::size_t N>
static_matrix<T, M, N> operator+(static_matrix<T, M, N> const& in) {
  return in;
}

template<typename T, std::size_t M, std::size_t N>
static_matrix<T, M, N> operator-(static_matrix<T, M, N> const& in) {
  static_matrix<T, M, N> out;
  for (typename static_matrix<T, M, N>::size_type i = 0; i < N; ++i) out[i] = -in[i];
  return out;
}

template<typename T, std::size_t M, std::size_t N>
auto operator+(static_matrix<T, M, N> lhs, static_matrix<T, M, N> const& rhs) {
  lhs += rhs;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N>
auto operator-(static_matrix<T, M, N> lhs, static_matrix<T, M, N> const& rhs) {
  lhs -= rhs;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N, typename Scalar>
auto operator*(static_matrix<T, M, N> lhs, Scalar const& fac) {
  lhs *= fac;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N, typename Scalar>
auto operator*(Scalar const& fac, static_matrix<T, M, N> lhs) {
  lhs *= fac;
  return lhs;
}

template<typename T, std::size_t M, std::size_t N, typename Scalar>
auto operator/(static_matrix<T, M, N> lhs, Scalar const& fac) {
  lhs /= fac;
  return lhs;
}

// Swap
template<typename T, std::size_t M, std::size_t N>
void swap(static_matrix<T, M, N>& lhs, static_matrix<T, M, N>& rhs) { lhs.swap(rhs); }

// Equality operators
template<typename T, std::size_t M, std::size_t N>
constexpr bool operator==(static_matrix<T, M, N> const& lhs,
                          static_matrix<T, M, N> const& rhs) {
  for (std::size_t i = 0; i < M; ++i)
    for (std::size_t j = 0; j < N; ++j)
      if (lhs(i, j) != rhs(i, j)) return false;
  return true;
}
template<typename T, std::size_t M, std::size_t N>
constexpr bool operator!=(static_matrix<T, M, N> const& lhs,
                          static_matrix<T, M, N> const& rhs) {
  return !(lhs == rhs);
}

} // namespace fivo

#endif // FIVO_MATRIX_H
