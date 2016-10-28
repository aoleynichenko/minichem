#include "QS_RuntimeError.h"
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

void QS_Scope::set(std::string name, QS_Object* obj)
{
  table[name] = obj;
}

QS_Object* QS_Scope::get(std::string name)
{
  auto it = table.find(name);
  if (it != table.end())
    return it->second;
  if (!top)
    throw QS_RuntimeError("undefined name: " + name);
  else
    return top->get(name);
}

  } // namespace qscript
} //namespace minichem
