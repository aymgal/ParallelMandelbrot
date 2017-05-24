#include <algorithm>

#include "grid.hh"

Grid::Grid(int m, int n) 
: m_nx(m), m_ny(n), m_storage(m * n) {
  clear();
}

void Grid::clear() {
  std::fill(m_storage.begin(), m_storage.end(), 0.0);
}

void Grid::resize(int nx, int ny) {
  m_nx = nx;
  m_ny = ny;
  m_storage.resize(nx * ny);
}

int Grid::nx() const { return m_nx; }
int Grid::ny() const { return m_ny; }

dfloat Grid::get_min() const {
  return *std::min_element(m_storage.begin(), m_storage.end());
}

dfloat Grid::get_max() const {
  return *std::max_element(m_storage.begin(), m_storage.end());
}
