#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>

#include "OutputStream.h"

namespace minichem {

using std::string;
using std::invalid_argument;

const string OutputStream::NEW    = "w";
const string OutputStream::APPEND = "a";
OutputStream OutputStream::STDOUT = OutputStream();
OutputStream OutputStream::STDERR = OutputStream(stderr);

OutputStream::OutputStream()
{
  file_ = stdout;
  isFlushingEnabled_ = false;
}

OutputStream::OutputStream(std::FILE* to)
{
  file_ = stderr;
  isFlushingEnabled_ = true;
}

OutputStream::OutputStream(string fileName, string mode)
{
  bindToFile(fileName, mode);
  isFlushingEnabled_ = false;
}

OutputStream* OutputStream::getStdout()
{
  return &STDOUT;
}

OutputStream* OutputStream::getStderr()
{
  return &STDERR;
}

void OutputStream::setFlushingEnabled(bool enabled)
{
  isFlushingEnabled_ = enabled;
  if (file_)
    fflush(file_);
}

void OutputStream::redirectToAnotherFile(std::string fileName, std::string mode)
{
  close();
  bindToFile(fileName, mode);
}

void OutputStream::close()
{
  if (file_)
    fclose(file_);
  file_ = 0;
}

// Object-oriented wrapper for fprintf()
// %b - print boolean
// others - see K&R 7.2
// C++ std::strings are not yet allowed
int OutputStream::printf(string format, ...)
{
  va_list ap;
  char* p, *sval;
  char *fmt = const_cast<char *>(format.c_str());
  int ival, count = 0;
  double dval;

  va_start(ap, fmt);
  for (p = fmt; *p; p++) {
    if (*p != '%') {
      fputc(*p, file_);
      continue;
    }
    char minfmt[16];
    char* mfp = minfmt;
    while (!isalpha(*p)) {
      *mfp++ = *p++;
    }
    *mfp++ = *p;
    *mfp = '\0';

    void *anyptr = 0;
    bool boolean;
    switch (*p) {
      case 'b':  // _b_ool
        boolean = va_arg(ap, int);
        fprintf(file_, "%s", (boolean == 1) ? "true" : "false");
        break;
      case 'd':
      case 'i':
      case 'o':
      case 'c':
        ival = va_arg(ap, int);
        fprintf(file_, minfmt, ival);
        break;
      case 'f':
      case 'e':
      case 'E':
      case 'g':
      case 'G':
        dval = va_arg(ap, double);
        fprintf(file_, minfmt, dval);
        break;
      case 's':
        sval = va_arg(ap, char *);
        fprintf(file_, minfmt, sval);
        break;
      case 'p':
        anyptr = va_arg(ap, void *);
        fprintf(file_, minfmt, anyptr);
      default:
        putchar(*p);
        break;
    }
    count++;
  }
  va_end(ap); /* очистка, когда все сделано */

  if (isFlushingEnabled_)
    fflush(file_);
  return count;
}

void OutputStream::println()
{
  this->printf("\n");  
}

OutputStream::~OutputStream()
{
  if (file_)
    fclose(file_);
}

void OutputStream::bindToFile(string fileName, string mode)
{
  if (mode == "a" || mode == "w") {
    file_ = fopen(fileName.c_str(), mode.c_str());
    if (!file_)
      throw invalid_argument("cannot open file " + fileName);
  }
  else
    throw invalid_argument("wrong mode: " + mode + " (only \"w\" or \"a\" are \
                            allowed)");
}


} // namespace minichem
