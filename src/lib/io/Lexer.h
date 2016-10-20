#ifndef _LEXER_H_INCLUDED
#define _LEXER_H_INCLUDED

#include "Token.h"

namespace minichem {

class Lexer
{
public:
	Lexer();

	Token get();
	Token match(Token t, int type);
	int   getint();
	double getdouble();
	void putback(Token t);
	void ignore(char c);
	void setInput(std::istream* input);
private:
	bool full_;
	Token buffer_;
	std::istream* inp_;
};

} // namespace minichem

#endif // _LEXER_H_INCLUDED
