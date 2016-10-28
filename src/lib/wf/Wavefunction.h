#ifndef _WAVEFUNCTION_H_INCLUDED
#define _WAVEFUNCTION_H_INCLUDED

#include <memory>

namespace minichem {

class Wavefunction {
public:
  virtual double energy() = 0;
  virtual ~Wavefunction();
protected:
  Wavefunction();
};

typedef std::shared_ptr<Wavefunction> Wf_ptr;

} // namespace minichem

#endif // _WAVEFUNCTION_H_INCLUDED
