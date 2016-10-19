#ifndef _TOKEN_H_INCLUDED
#define _TOKEN_H_INCLUDED

#include <ostream>
#include <string>

namespace minichem {

class Token
{
public:
  static const int TT_EOF;
  static const int TT_NOTHING;
  static const int TT_NUMBER;
  static const int TT_WORD;

  static const int TT_KW_FUN;
  static const int TT_KW_MOL;
  static const int TT_KW_VAR;

  Token(int toktype);
  Token(int toktype, double val);
  Token(int toktype, std::string str);

  std::string toString() const;

  int ttype;
  std::string sval;
  double dval;
private:
  int lookupKeyword(std::string word);
};

std::ostream& operator<<(std::ostream& out, const Token& t);

} // namespace minichem

#endif // _TOKEN_H_INCLUDED
