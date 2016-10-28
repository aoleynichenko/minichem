/* Represents scope of view with names.
 * In fact, QS_Scope is a sequence containers map<name,QS_Object>, each
 * represents one scope of view (inside function, loop, for example).
 *
 * QuantumScript Engine contains only one instance of QS_Scope object
 * (top-level scope).
 * Each scope has a pointer to a wider scope. Hence, we should use recursive
 * algorithm to get an object by its name:
 *
 *  GetObject(name):
 *   if name in table.keys():
 *     return table[name]
 *   if top == NULL:
 *     error 'not found'
 *   return top->GetObject(name)  // recursion
 */

#ifndef _QS_SCOPE_H_INCLUDED
#define _QS_SCOPE_H_INCLUDED

#include <map>
#include <string>

#include "QS_Object.h"

namespace minichem {
  namespace qscript {

class QS_Scope {
public:
  QS_Scope();
  ~QS_Scope();
  void newScope();
  QS_Object* operator[](std::string name);
private:
  QS_Scope* top;  // scope 'this' is nested into 'top'
  std::map<std::string, QS_Object> table;
};

  }
} //namespace minichem

#endif // _QS_SCOPE_H_INCLUDED
