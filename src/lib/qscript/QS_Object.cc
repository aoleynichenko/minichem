#include <sstream>
#include <string>

#include "QS_Object.h"

namespace minichem {
  namespace qscript {

using std::ostringstream;
using std::string;

// protected contructor
QS_Object::QS_Object()
{
}

// by default, returns address of this Object
std::string QS_Object::toString() const
{
  const void* address = static_cast<const void*>(this);
  ostringstream ss;
  ss << address;
  return ss.str();
}

QS_Object::~QS_Object()
{

}

  } // namespace qscipt
} // namespace minichem
