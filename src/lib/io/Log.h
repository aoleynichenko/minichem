#ifndef _LOG_H_INCLUDED
#define _LOG_H_INCLUDED

#include <cstdio>
#include <string>

namespace minichem {

class Log {
public:
	Log(std::string logname);
	void log(std::string format, ...);
	~Log();
private:
	FILE* logfile_m;
	std::string currentDateToString();
};

} // namespace minichem

#endif // _LOG_H_INCLUDED