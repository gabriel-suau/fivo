#ifndef FIVO_TRAITS_H
#define FIVO_TRAITS_H

#include <type_traits>

namespace fivo {

namespace traits {

template<typename Derived, typename Base>
struct is_derived {
  static constexpr bool value =
    std::is_base_of_v<Base, Derived> &&
    std::is_convertible_v<const volatile Derived*, const volatile Base*>;
};

} // namespace traits

} // namespace fivo

#endif // FIVO_TRAITS_H
