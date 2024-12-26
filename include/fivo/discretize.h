#ifndef FIVO_DISCRETIZE_H
#define FIVO_DISCRETIZE_H

#include <fivo/state.h>
#include <fivo/mesh.h>

#include <utility>

namespace fivo {

template<typename Func>
auto discretize(Mesh const& mesh, Func&& func) {
  using state_type = decltype(func(std::declval<double>()));
  global_state<state_type> X(mesh.nx());
  for (int i = 0; i < mesh.nx(); ++i) { X[i] = func(mesh.cell_center(i)); }
  return X;
}

} // namespace fivo

#endif // FIVO_DISCRETIZE_H
