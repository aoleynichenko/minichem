#ifndef _QS_OBJECT_H_INCLUDED
#define _QS_OBJECT_H_INCLUDED

#include <string>

namespace minichem {
  namespace qscript {

class QS_Object {
public:
  virtual std::string toString() const;
  virtual ~QS_Object();
protected:
  QS_Object();
};

  } // namespace qscipt
} // namespace minichem

#endif // _QS_OBJECT_H_INCLUDED
