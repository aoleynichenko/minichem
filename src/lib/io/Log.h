#ifndef _LOG_H_INCLUDED
#define _LOG_H_INCLUDED

#include <cstdio>
#include <string>

namespace minichem {

std::string detectCmdLog(int argc, char *argv[]);


class Log {
public:

	Log(std::string logname);
	void log(std::string format, ...);
	~Log();

	inline std::string getLogName()
	{
		return logName;
	}
private:
	std::string logName;
	FILE* logfile_m;
	std::string currentDateToString();
};

} // namespace minichem

#endif // _LOG_H_INCLUDED
