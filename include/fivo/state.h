#ifndef FIVO_STATE_H
#define FIVO_STATE_H

#include <fivo/array.h>
#include <fivo/vector.h>

namespace fivo {

/* LOCAL STATE */
template<typename T, std::size_t NUM>
using state = array<T, NUM>;

/* GLOBAL STATE */
template<typename State>
using global_state = vector<State>;

} // namespace fivo

#endif // FIVO_STATE_H
