#ifndef FIVO_ARRAY_H
#define FIVO_ARRAY_H

#include <type_traits>
#include <tuple>

namespace fivo {

/** ARRAY HELPER */
template<typename T, std::size_t N>
struct array {
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
  using reverse_iterator       = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  // Data
  T storage[N] = {};

  // Element access
  constexpr const_reference at(size_type i) const {
    if (i >= size())
      throw std::out_of_range("fivo::array: Index out of range.");
    return storage[i];
  }
  constexpr reference at(size_type i) {
    if (i >= size())
      throw std::out_of_range("fivo::array: Index out of range.");
    return storage[i];
  }
  constexpr const_reference operator[](size_type i) const noexcept { return storage[i]; }
  constexpr       reference operator[](size_type i)       noexcept { return storage[i]; }
  constexpr const_reference front() const noexcept { return storage[0]; }
  constexpr       reference front()       noexcept { return storage[0]; }
  constexpr const_reference back()  const noexcept { return storage[N - 1]; }
  constexpr       reference back()        noexcept { return storage[N - 1]; }
  constexpr const_pointer   data()  const noexcept { return &storage[0]; }
  constexpr       pointer   data()        noexcept { return &storage[0]; }

  // Iterators
  constexpr       iterator begin()        noexcept { return       iterator(data()); }
  constexpr const_iterator begin()  const noexcept { return const_iterator(data()); }
  constexpr const_iterator cbegin() const noexcept { return const_iterator(data()); }
  constexpr       iterator end()          noexcept { return       iterator(data() + N); }
  constexpr const_iterator end()    const noexcept { return const_iterator(data() + N); }
  constexpr const_iterator cend()   const noexcept { return const_iterator(data() + N); }
  constexpr       reverse_iterator rbegin()        noexcept { return       reverse_iterator(end()); }
  constexpr const_reverse_iterator rbegin()  const noexcept { return const_reverse_iterator(end()); }
  constexpr const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(end()); }
  constexpr       reverse_iterator rend()          noexcept { return       reverse_iterator(begin() + N); }
  constexpr const_reverse_iterator rend()    const noexcept { return const_reverse_iterator(begin() + N); }
  constexpr const_reverse_iterator crend()   const noexcept { return const_reverse_iterator(begin() + N); }

  // Capacity
  constexpr bool empty() const noexcept { return size() == 0; }
  constexpr size_type size() const noexcept { return N; }
  constexpr size_type max_size() const noexcept { return N; }

  // Other
  constexpr void fill(T const& value) { for (size_type i = 0; i < N; ++i) storage[i] = value; }
  constexpr void swap(array& other) noexcept(std::is_nothrow_swappable_v<T>) {
    array tmp = other;
    other = *this;
    *this = tmp;
  }

  // Arithmetic operators
  constexpr inline auto& operator+=(array<T, N> const& rhs) {
    for (size_type i = 0; i < N; ++i) storage[i] += rhs[i];
    return *this;
  }
  constexpr inline auto& operator-=(array<T, N> const& rhs) {
    for (size_type i = 0; i < N; ++i) storage[i] -= rhs[i];
    return *this;
  }
  template<typename Scalar,
           std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
  constexpr inline auto& operator*=(Scalar const& factor) {
    for (size_type i = 0; i < N; ++i) storage[i] *= factor;
    return *this;
  }
  template<typename Scalar,
           std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
  constexpr inline auto& operator/=(Scalar const& factor) {
    for (size_type i = 0; i < N; ++i) storage[i] /= factor;
    return *this;
  }
};

template<typename T, std::size_t N>
constexpr inline array<T, N> operator+(array<T, N> const& in) {
  return in;
}

template<typename T, std::size_t N>
constexpr inline array<T, N> operator-(array<T, N> const& in) {
  array<T, N> out;
  for (typename array<T, N>::size_type i = 0; i < N; ++i) out[i] = -in[i];
  return out;
}

template<typename T, std::size_t N>
constexpr inline array<T, N> operator+(array<T, N> lhs, array<T, N> const& rhs) {
  lhs += rhs;
  return lhs;
}

template<typename T, std::size_t N>
constexpr inline array<T, N> operator-(array<T, N> lhs, array<T, N> const& rhs) {
  lhs -= rhs;
  return lhs;
}

template<typename T, std::size_t N, typename Scalar,
         std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
constexpr inline array<T, N> operator*(array<T, N> lhs, Scalar const& fac) {
  lhs *= fac;
  return lhs;
}

template<typename T, std::size_t N, typename Scalar,
         std::enable_if_t<std::is_arithmetic<Scalar>::value, bool> = false>
constexpr inline array<T, N> operator*(Scalar const& fac, array<T, N> lhs) {
  lhs *= fac;
  return lhs;
}

template<typename T, std::size_t N, typename Scalar>
constexpr inline array<T, N> operator/(array<T, N> lhs, Scalar const& fac) {
  lhs /= fac;
  return lhs;
}

// CTAD
template<typename T, typename... U>
array(T, U...) -> array<std::enable_if_t<(std::is_same_v<T, U> && ...)>,
                        1 + sizeof...(U)>;

// Equality operators
template<typename T, std::size_t N>
constexpr inline bool operator==(array<T, N> const& lhs, array<T, N> const& rhs) {
  for (std::size_t i = 0; i < N; ++i) { if (lhs[i] != rhs[i]) return false; }
  return true;
}
template<typename T, std::size_t N>
constexpr inline bool operator!=(array<T, N> const& lhs, array<T, N> const& rhs) {
  return !(lhs == rhs);
}
template<typename T, std::size_t N>
constexpr inline bool operator<(array<T, N> const& lhs, array<T, N> const& rhs) {
  return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}
template<typename T, std::size_t N>
constexpr inline bool operator>(array<T, N> const& lhs, array<T, N> const& rhs) {
  return rhs < lhs;
}
template<typename T, std::size_t N>
constexpr inline bool operator<=(array<T, N> const& lhs, array<T, N> const& rhs) {
  return !(lhs > rhs);
}
template<typename T, std::size_t N>
constexpr inline bool operator>=(array<T, N> const& lhs, array<T, N> const& rhs) {
  return !(lhs < rhs);
}

} // namespace fivo

namespace std {

template<typename T, std::size_t N>
void swap(::fivo::array<T, N>& lhs, ::fivo::array<T, N>& rhs) noexcept(noexcept(lhs.swap(rhs))) {
  lhs.swap(rhs);
}

} // namespace std

// Partial specialization of std::tuple_size and std::tuple_element for
// array to structured bindings
template<typename T, std::size_t N>
struct std::tuple_size<fivo::array<T, N>>
  : std::integral_constant<std::size_t, N> {};

template<std::size_t I, typename T, std::size_t N>
struct std::tuple_element<I, fivo::array<T, N>> {
  using type = T;
};

namespace fivo {

template<std::size_t I, class T, std::size_t N>
constexpr T& get(array<T, N>& a) noexcept {
  return a[I];
}
template<std::size_t I, class T, std::size_t N>
constexpr T const& get(array<T, N> const& a) noexcept {
  return a[I];
}
template<std::size_t I, class T, std::size_t N>
constexpr T&& get(array<T, N>&& a) noexcept {
  return std::move(a[I]);
}
template<std::size_t I, class T, std::size_t N>
constexpr T const&& get(array<T, N> const&& a) noexcept {
  return std::move(a[I]);
}

} // namespace fivo

#endif // FIVO_ARRAY_H
