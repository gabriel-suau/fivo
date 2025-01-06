#ifndef FIVO_STATE_H
#define FIVO_STATE_H

#include <fivo/array.h>
#include <fivo/vector.h>
#include <fivo/matrix.h>

namespace fivo {

/* LOCAL STATE */
template<typename T, std::size_t NUM>
using state = array<T, NUM>;

/* GEOMETRIC VECTOR */
template<typename T, std::size_t DIM>
using vec = array<T, DIM>;

/* LOCAL FLUX */
template<typename State, std::size_t DIM>
using flux = static_matrix<typename State::value_type,
                           std::tuple_size<State>::value, DIM>;

/* GLOBAL STATE */
template<typename State>
using global_state = vector<State>;

/* GLOBAL FLUX */
template<typename Flux>
using global_flux = vector<Flux>;

} // namespace fivo

#endif // FIVO_STATE_H
