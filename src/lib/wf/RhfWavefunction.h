#ifndef _RHF_WAVEFUNCTION_H_INCLUDED
#define _RHF_WAVEFUNCTION_H_INCLUDED

#include "Wavefunction.h"

namespace minichem {

class RhfWavefunction {
public:
  RhfWavefunction();
  double energy();
  ~RhfWavefunction();
};

} // namespace minichem

#endif // _RHF_WAVEFUNCTION_H_INCLUDED
