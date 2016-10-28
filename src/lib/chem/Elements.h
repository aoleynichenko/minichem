#ifndef _ELEMENTS_H_INCLUDED
#define _ELEMENTS_H_INCLUDED

#include <string>

namespace minichem {

struct Element {
  int Z;
  double mass;
  std::string sym;
  std::string name;
};

// interface to static std::vector<Element> info
struct Elements {
  static Element& getElementByZ(int z);
  static Element& getElementBySym(std::string sym);
  static int sym2charge(std::string sym);
  static std::string charge2sym(int charge);
};

} // namespace minichem

#endif  // _ELEMENTS_H_INCLUDED
