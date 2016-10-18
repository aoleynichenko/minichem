#include <cctype>
#include <iostream>
#include <string>

#include "Lexer.h"
#include "Token.h"

namespace minichem {

using std::cin;
using std::cout;
using std::endl;
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
