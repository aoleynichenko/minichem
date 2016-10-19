#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "Lexer.h"
#include "Token.h"
#include "../except/SyntaxError.h"

namespace minichem {

using std::cin;
using std::cout;
using std::endl;
using std::ostringstream;
using std::string;

Lexer::Lexer()
  : full_(false), buffer_('\0'), inp_(&cin)
{
}

Token Lexer::get()
{
  if (full_) {
    full_ = false;
    return buffer_;
  }

  char ch;
  *inp_ >> ch;
  if (!*inp_)
    return Token(Token::TT_EOF);
  if (isdigit(ch)) {
    inp_->putback(ch);
    double val;
    *inp_ >> val;
    return Token(Token::TT_NUMBER, val);
  }
  else if (isalpha(ch)) {
    string word;
    while (isalnum(ch) || ch == '_') {
      word += ch;
      ch = inp_->get();
    }
    inp_->putback(ch);
    return Token(Token::TT_WORD, word);
  }
  else
    return Token(ch);
}

Token Lexer::match(Token t, int type)
{
  if (t.ttype != type) {
    ostringstream errmsg;
    errmsg << "expected token type: {int: " << type << ", char: " << (char) type << "}, "
      << "but found " << t.toString();
    throw SyntaxError(errmsg.str());
  }
  return get();
}

int Lexer::getint()
{
  Token t = get();
  if (t.ttype != Token::TT_NUMBER)
    throw SyntaxError("expected integer number, but found " + t.toString());
  // is t.dval integer?
  if (fabs(t.dval) != fabs((int) t.dval))
    throw SyntaxError("expected integer number, but found double");
  return (int) t.dval;
}

void Lexer::putback(Token t)
{
  full_ = true;
  buffer_ = t;
}

void Lexer::ignore(char c)
{

}

void Lexer::setInput(std::istream* input)
{
  inp_ = input;
  full_ = false;
}

} // namespace minichem
