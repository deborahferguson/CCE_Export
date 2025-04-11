#include "utils.hh"

namespace CCE_export {

int l_m_to_index(int l, int m) { 
  return l * l + l + m; 
}

} //namespace CCE_export
