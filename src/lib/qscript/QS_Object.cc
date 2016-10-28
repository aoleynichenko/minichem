#include <sstream>
#include <string>
#include <vector>

#include "QS_Object.h"

namespace minichem {
  namespace qscript {

using std::ostringstream;
using std::string;
using std::vector;

QS_Object::QS_Object()
  : type(TYPE_OBJ)
{
}

// by default, returns string in format {key1: value1, key2: value2, ... }
// if object has FUNction str(), invokes this function (str === toString)
string QS_Object::toString() const
{
  ostringstream ss;
  vector<string> strs;
  ss << "{";
  for (const auto& kv : fields_)
    strs.push_back(kv.first + ": " + kv.second->toString());
  for (vector<string>::size_type i = 0; i < strs.size(); i++)
    ss << strs[i] << ((i < strs.size() - 1) ? ", " : "");
  ss << "}";
  return ss.str();
}

string QS_Object::getTypeString() const
{
  return "object";
}

QS_Object::~QS_Object()
{

}

  } // namespace qscipt
} // namespace minichem
