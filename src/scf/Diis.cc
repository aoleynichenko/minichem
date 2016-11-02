#include "Diis.h"

namespace minichem {

Diis::Diis(int diisbas, const Matrix& S)
  : bas_m(diisbas), S_m(S)
{
}

Diis::~Diis()
{
}

} //namespace minichem
