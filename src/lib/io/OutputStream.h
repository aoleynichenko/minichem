#ifndef _OUTPUT_STREAM_H_INCLUDED
#define _OUTPUT_STREAM_H_INCLUDED

#include <cstdio>
#include <string>

namespace minichem {

class OutputStream
{
public:
  static const std::string NEW;
  static const std::string APPEND;

  OutputStream(std::string fileName, std::string mode = NEW);

  static OutputStream* getStdout();
  static OutputStream* getStderr();

  void setFlushingEnabled(bool enabled);
  void redirectToAnotherFile(std::string fileName, std::string mode = NEW);
  void close();
  int  printf(std::string format, ...);
  void println();

  ~OutputStream();
protected:
  OutputStream();
  OutputStream(std::FILE* to);
private:
  bool isFlushingEnabled_;
  std::FILE *file_;
  static OutputStream STDOUT;
  static OutputStream STDERR;

  void bindToFile(std::string fileName, std::string mode);
};

} // namespace minichem

#endif // _OUTPUT_STREAM_H_INCLUDED