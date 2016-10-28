#include <ostream>
#include <sstream>
#include <string>

#include "Token.h"

namespace minichem {

using std::ostream;
using std::ostringstream;
using std::string;

const int Token::TT_EOF     = -1;
const int Token::TT_EOL     = -2;
const int Token::TT_NOTHING = -3;
const int Token::TT_NUMBER  = -4;
const int Token::TT_WORD    = -5;
const int Token::TT_QUOTE    = -6;

const int Token::TT_KW_FUN   = -10;
const int Token::TT_KW_MOL   = -11;
const int Token::TT_KW_BASIS = -12;
const int Token::TT_KW_VAR   = -13;
const int Token::TT_KW_TASK  = -14;

Token::Token(int toktype)
  : ttype(toktype), sval(""), dval(0.0)
{
}

Token::Token(int toktype, double val)
  : ttype(toktype), sval(""), dval(val)
{
}

Token::Token(int toktype, string str)
  : sval(str), dval(0.0)
{
  ttype = lookupKeyword(str);
}

string Token::toString() const
{
  ostringstream s;
  s << "[";
  if (ttype == TT_EOF)
    s << "END OF FILE]";
  else if (ttype == TT_EOL)
    s << "END OF LINE]";
  else if (ttype == TT_NUMBER)
    s << "NUMBER  | " << dval << "]";
  else if (ttype == TT_WORD)
    s << "WORD    | " << sval << "]";
  else if (ttype <= TT_KW_FUN)
    s << "KEYWORD | " << sval << "]";
  else
    s << (char) ttype << "]";
  return s.str();
}

int Token::lookupKeyword(string word)
{
  if (word == "var")
    return TT_KW_VAR;
  else if (word == "mol")
    return TT_KW_MOL;
  else if (word == "basis")
    return TT_KW_BASIS;
  else if (word == "fun")
    return TT_KW_FUN;
  else if (word == "task")
    return TT_KW_TASK;
  else
    return TT_WORD;
}

ostream& operator<<(ostream& out, const Token& t)
{
  out << t.toString();
  return out;
}

} // namespace minichem
