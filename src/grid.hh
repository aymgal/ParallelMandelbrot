#ifndef GRID_HH
#define GRID_HH

#include <vector>

#include "gvars.hh"

class Grid {
public:
  Grid(int m, int n);

  /// access the value [i][j] of the grid
  inline dfloat & operator()(int ix, int iy) {
    return m_storage[ix * m_ny + iy];
  }
  inline const dfloat & operator()(int ix, int iy) const {
    return m_storage[ix * m_ny + iy];
  }

  void clear();

  void resize(int nx, int ny);

  int nx() const;
  int ny() const;

  dfloat get_min() const;
  dfloat get_max() const;

private:
  int m_nx, m_ny;
  std::vector<dfloat> m_storage;
};

#endif /* GRID_HH */
