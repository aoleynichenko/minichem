/* Note: this class should be independent on LIBINT to provide compatibility
 * with minichem's native AO integrals engine.
 */

#ifndef _BASIS_SET_H_INCLUDED
#define _BASIS_SET_H_INCLUDED

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace minichem {

class BasisSet {
public:
  // LBlock struct designed to be as close as possible to libint2::Shell to
  // simplify generation of atom-centered basis set (== vector<Shell>)
  struct LBlock {
    int l_;
    bool cart_;
    std::vector<double> alpha_;
    std::vector<std::vector<double>> contr_; // contractions
  };

  typedef std::vector<LBlock> BasisTemplate_t; // basis set template for one element
  static std::string am2string(int L);  // L -> {S, P, D, F...}
  static int parseAngmom(std::string am);  // {S, P, D, F... -> L}

  void addLBlock(std::string elemSym, LBlock block);
  std::string toString() const;
private:
  std::map<int, BasisTemplate_t> set_;
};

typedef std::shared_ptr<BasisSet> BasisSet_ptr;

} // namespace minichem

#endif // _BASIS_SET_H_INCLUDED
