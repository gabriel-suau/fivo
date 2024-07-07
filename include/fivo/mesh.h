#ifndef FIVO_MESH_H
#define FIVO_MESH_H

#include <vector>

namespace fivo {

struct Mesh {
  Mesh() = default;
  Mesh(double xmin, double xmax, int nx)
    : m_xmin(xmin), m_xmax(xmax), m_dx((xmax - xmin) / nx), m_nx(nx)
  {}

  auto xmin() const { return m_xmin; }
  auto xmax() const { return m_xmax; }
  auto dx() const { return m_dx; }
  auto nx() const { return m_nx; }

  auto cell_center(int i) const { return xmin() + (i + 0.5) * dx(); }

  void xmin(double value) {
    m_xmin = value;
    m_dx = (m_xmax - m_xmin) / m_nx;
  }
  void xmax(double value) {
    m_xmax = value;
    m_dx = (m_xmax - m_xmin) / m_nx;
  }
  void nx(int value) {
    m_nx = value;
    m_dx = (m_xmax - m_xmin) / m_nx;
  }

private:
  double m_xmin, m_xmax, m_dx;
  int m_nx;
};

} // namespace fivo

#endif // FIVO_MESH_SH
