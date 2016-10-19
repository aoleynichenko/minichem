#ifndef _SYNTAX_ERROR_H_INCLUDED
#define _SYNTAX_ERROR_H_INCLUDED

#include <stdexcept>
#include <string>

namespace minichem {

class SyntaxError : public std::runtime_error
{
public:
  SyntaxError(const std::string& message)
    : std::runtime_error(message) {}
};

} // namespace minichem

#endif // _SYNTAX_ERROR_H_INCLUDED
