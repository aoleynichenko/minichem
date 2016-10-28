// in QuantumScript, all is object

#ifndef _QS_OBJECT_H_INCLUDED
#define _QS_OBJECT_H_INCLUDED

#include <map>
#include <string>

namespace minichem {
  namespace qscript {

class QS_Object {
public:
  enum Type {TYPE_OBJ, TYPE_STR, TYPE_FUN, TYPE_INT, TYPE_DBL, TYPE_ARR,
             TYPE_MOL, TYPE_BAS, TYPE_WFN};
  Type type;

  QS_Object();
  virtual std::string toString() const;
  virtual std::string getTypeString() const;
  virtual ~QS_Object();
private:
  std::map<std::string, QS_Object*> fields_;
};

  } // namespace qscipt
} // namespace minichem

#endif // _QS_OBJECT_H_INCLUDED
