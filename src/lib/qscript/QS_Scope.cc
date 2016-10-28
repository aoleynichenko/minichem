#include "QS_Scope.h"

namespace minichem {
  namespace qscript {

QS_Scope::QS_Scope()
{
  top = nullptr;
}

QS_Scope::~QS_Scope()
{
}

void QS_Scope::newScope()
{

}

QS_Object* QS_Scope::operator[](std::string name)
{
  
}

  } // namespace qscript
} //namespace minichem
