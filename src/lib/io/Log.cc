#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <string>

#include "Log.h"

using std::invalid_argument;
using std::string;

namespace minichem {

Log::Log(string logname)
{
	logfile_m = fopen(logname.c_str(), "a");
	if (!logfile_m)
		throw invalid_argument("Cannot create log file " + logname);
	log("Log created");
}

// this method is just revised OutputStream::printf(string format, ...)
// automatically puts '\n' at the end of line
void Log::log(string format, ...)
{
	va_list ap;
	char* p, *sval;
	char *fmt = const_cast<char *>(format.c_str());
	int ival;
 	double dval;

 	// start with date:
 	fprintf(logfile_m, "[%s] ", currentDateToString().c_str());
 	// loop over format string
 	va_start(ap, fmt);
 	for (p = fmt; *p; p++) {
		if (*p != '%') {
			fputc(*p, logfile_m);
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
			fprintf(logfile_m, "%s", (boolean == 1) ? "true" : "false");
			break;
		case 'd':
		case 'i':
		case 'o':
		case 'c':
			ival = va_arg(ap, int);
			fprintf(logfile_m, minfmt, ival);
			break;
		case 'f':
		case 'e':
		case 'E':
		case 'g':
		case 'G':
			dval = va_arg(ap, double);
			fprintf(logfile_m, minfmt, dval);
			break;
		case 's':
			sval = va_arg(ap, char *);
			fprintf(logfile_m, minfmt, sval);
			break;
		case 'p':
			anyptr = va_arg(ap, void *);
			fprintf(logfile_m, minfmt, anyptr);
		default:
			putchar(*p);
			break;
		}
	}
	va_end(ap);

	fputc('\n', logfile_m);
	fflush(logfile_m);
}

Log::~Log()
{
	fclose(logfile_m);
}

// source: http://stackoverflow.com/questions/27658477/date-in-log-file-using-ofstream
string Log::currentDateToString()
{  
    time_t rawtime;
    struct tm* timeInfo;
    char buffer[80];

    time(&rawtime);
    timeInfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeInfo);
    string dateString(buffer);

    return dateString;
}

} // namespace minichem
