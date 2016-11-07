#ifndef _UHF_WAVEFUNCTION_H_INCLUDED
#define _UHF_WAVEFUNCTION_H_INCLUDED

#include "Wavefunction.h"

namespace minichem {

class UhfWavefunction {
public:
  UhfWavefunction();
  double energy();
  ~UhfWavefunction();
};

} // namespace minichem

#endif // _UHF_WAVEFUNCTION_H_INCLUDED
