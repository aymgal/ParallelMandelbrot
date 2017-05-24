#ifndef BUFFER_HH
#define BUFFER_HH

#include <memory>

#include "grid.hh"

class Buffer {
public:
  Buffer(int m, int n);

  Grid & storage();
  void resize(int m, int n);

private:
  std::unique_ptr<Grid> m_pstorage;
};

#endif /* BUFFER_HH */
