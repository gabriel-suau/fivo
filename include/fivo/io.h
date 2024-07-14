#ifndef FIVO_IO_H
#define FIVO_IO_H

#include "state.h"
#include "mesh.h"
#include <string>
#include <fstream>
#include <tuple>
#include <ctime>

namespace fivo {

namespace impl {

template<typename Func, typename Tuple, std::size_t... Is>
void apply_impl(Func const& func, Tuple const& t, std::index_sequence<Is...>) {
  (func(std::get<Is>(t)), ...);
}

template<typename Func, typename... Args>
void apply(Func const& func, std::tuple<Args...> const& t) {
  apply_impl(func, t, std::index_sequence_for<Args...>{});
}

} // namespace impl

class IOManager {
public:
  IOManager(std::string const& basename, int save_frequency, Mesh const& mesh)
    : m_basename(basename), m_save_frequency(save_frequency), m_mesh(mesh)
  {}

  auto const& basename() const { return m_basename; } 
  auto save_frequency() const { return m_save_frequency; }
  auto const& mesh() const { return m_mesh; }
  void basename(std::string const& new_basename) { m_basename = new_basename; }
  void save_frequency(int const& new_frequency) { m_save_frequency = new_frequency; }
  void mesh(Mesh const& new_mesh) { m_mesh = new_mesh; }

  template<typename State, typename... Quantities>
  void save_state(typename State::value_type const& t, int niter,
                  global_state<State> const& X,
                  std::tuple<Quantities...> const& quantities) const {
    std::string filename = m_basename + "_" + std::to_string(niter) + ".dat";
    std::ofstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Could not open output file " + filename);

    // Save header as comment
    auto const date_time = std::time(nullptr);
    char dtstr[100];
    std::strftime(dtstr, sizeof(dtstr), "%c", std::localtime(&date_time));
    file << "# fivo output\n"
         << "# " << dtstr << "\n"
         << "# at time t = " << t << "\n";

    // Save the quantities names
    file << "# x";
    impl::apply([&] (auto const& q) { file << " " << q.first; }, quantities);
    file << "\n";

    // Save the quantities values
    for (int i = 0; i < m_mesh.nx(); ++i) {
      auto const x = m_mesh.cell_center(i);
      file << x;
      impl::apply([&] (auto const& q) { file << " " << q.second(x, X[i]); }, quantities);
      file << "\n";
    }
  }

  template<typename State, typename... Quantities>
  void save_state(typename State::value_type const& t, int niter,
                  global_state<State> const& X,
                  Quantities... quantities) const {
    return save_state(t, niter, X, std::make_tuple(quantities...));
  }

private:
  std::string m_basename;
  int m_save_frequency;
  Mesh m_mesh;
};


} // namespace fivo

#endif // FIVO_IO_H
