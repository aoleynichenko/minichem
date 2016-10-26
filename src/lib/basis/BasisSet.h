/* Note: this class should be independent on LIBINT to provide compatibility
 * with minichem's native AO integrals engine.
 */

#ifndef _BASIS_SET_H_INCLUDED
#define _BASIS_SET_H_INCLUDED

#include <map>
#include <string>
#include <vector>

namespace minichem {

class BasisSet {
public:
  // LBlock struct designed to be as close as possible to libint2::Shell to
  // simplify generation of atom-centered basis set (== vector<Shell>)
  struct LBlock {
    int l_;
    std::vector<double> alpha_;
    std::vector<std::vector<double>> contr_; // contractions
  };

  typedef std::vector<LBlock> BasisTemplate_t; // basis set template for one element

  void addLBlock(std::string elemSym, LBlock block);
private:
  std::map<std::string, BasisTemplate_t> set_;
};

} // namespace minichem

#endif // _BASIS_SET_H_INCLUDED
