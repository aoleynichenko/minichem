#ifndef _QS_RUNTIME_ERROR_H_INCLUDED
#define _QS_RUNTIME_ERROR_H_INCLUDED

#include <stdexcept>
#include <string>

#include "QS_Object.h"

namespace minichem {
  namespace qscript {

class QS_RuntimeError : public std::runtime_error, QS_Object {
public:
  QS_RuntimeError(std::string);
};

  } // namespace qscript
} // namespace minichem

#endif // _QS_RUNTIME_ERROR_H_INCLUDED
