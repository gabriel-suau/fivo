#ifndef FIVO_TRAITS_H
#define FIVO_TRAITS_H

#include <type_traits>
#include <tuple>
#include <utility>

namespace fivo {

namespace traits {

template<typename Derived, typename Base>
struct is_derived {
  static constexpr bool value =
    std::is_base_of_v<Base, Derived> &&
    std::is_convertible_v<const volatile Derived*, const volatile Base*>;
};
template<typename Derived, typename Base>
constexpr bool is_derived_v = is_derived<Derived, Base>::value;

namespace impl {

template<typename Tuple, typename Func, std::size_t... Is>
void for_each_impl(Tuple&& tuple, Func&& func, std::index_sequence<Is...>) {
  (func(std::get<Is>(tuple)), ...);
}

} // namespace impl

template<typename Tuple, typename Func>
void for_each(Tuple&& tuple, Func&& func) {
  constexpr auto size = std::tuple_size<std::decay_t<Tuple>>::value;
  impl::for_each_impl(std::forward<Tuple>(tuple), std::forward<Func>(func),
                      std::make_index_sequence<size>());
}

// Helper to create a tuple with pairs (Cartesian product of elements)
namespace impl {

template<typename T, typename Tuple, std::size_t... Is>
auto create_pairs(T const& elem, Tuple const& t, std::index_sequence<Is...>) {
  return std::make_tuple(std::make_pair(elem, std::get<Is>(t))...);
}

template<typename Tuple, std::size_t... Is>
auto flatten_tuple_of_tuples_impl(Tuple const& t, std::index_sequence<Is...>) {
  return std::tuple_cat(std::get<Is>(t)...);
}

template<typename Tuple>
auto flatten_tuple_of_tuples(Tuple const& t) {
  constexpr auto N = std::tuple_size<Tuple>::value;
  return flatten_tuple_of_tuples_impl(t, std::make_index_sequence<N>{});
}

// Main Cartesian product function
template<typename... Ts1, typename... Ts2, std::size_t... Is>
auto cartesian_product_impl(std::tuple<Ts1...> const& t1,
                            std::tuple<Ts2...> const& t2,
                            std::index_sequence<Is...>) {
  return flatten_tuple_of_tuples(std::make_tuple
                                 (create_pairs
                                  (std::get<Is>(t1), t2,
                                   std::index_sequence_for<Ts2...>{})...));
}

} // namespace impl

template<typename... Ts1, typename... Ts2>
auto cartesian_product(std::tuple<Ts1...> const& t1, std::tuple<Ts2...> const& t2) {
  return impl::cartesian_product_impl(t1, t2, std::index_sequence_for<Ts1...>{});
}

} // namespace traits

} // namespace fivo

#endif // FIVO_TRAITS_H
