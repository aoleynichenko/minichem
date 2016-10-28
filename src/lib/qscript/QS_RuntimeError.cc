#include <stdexcept>
#include <string>

#include "QS_Object.h"
#include "QS_RuntimeError.h"

namespace minichem {
  namespace qscript {

using std::runtime_error;

QS_RuntimeError::QS_RuntimeError(std::string msg)
  : runtime_error(msg), QS_Object()
{
}

  } // namespace qscript
} // namespace minichem
