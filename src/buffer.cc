#include "buffer.hh"

Buffer::Buffer(int m, int n)
    : m_pstorage(new Grid(m, n)) {}

Grid & Buffer::storage() {
  return *m_pstorage;
}

void Buffer::resize(int m, int n) {
  m_pstorage->resize(m, n);
}