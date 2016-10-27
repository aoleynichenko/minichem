#ifndef _WAVEFUNCTION_H_INCLUDED
#define _WAVEFUNCTION_H_INCLUDED

namespace minichem {

class Wavefunction {
public:
  virtual double energy() = 0;
  virtual ~Wavefunction();
protected:
  Wavefunction();
};

} // namespace minichem

#endif // _WAVEFUNCTION_H_INCLUDED
