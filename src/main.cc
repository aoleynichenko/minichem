/* Minichem - open-source simple quantum chemistry code,
 * written in educational purposes.
 * 
 * Alexander Oleynichenko, Moscow State University, 2016
 * alexvoleynichenko@gmail.com
 *
 *  main.cc
 *  Kernel of minichem application.
 *  Invokes all another modules: input, hf, mp2, cc, geom, etc.
 */

#include <chrono>
#include <ctime>
#include <stdexcept>
#include <string>

#include "Kernel.h"
#include "minichem.h"
#include "lib/io/OutputStream.h"

using std::exception;
using std::string;

namespace minichem {

// GLOBAL DECLARATIONS
Log* mainlog = NULL;
std::chrono::high_resolution_clock::time_point startTimePoint;

} // namespace minichem

using namespace minichem;

void outputHeader()
{
	time_t t = time(0);
	OutputStream* out = OutputStream::getStdout();

	out->println();
	out->printf("                           **********************\n");
	out->printf("                           *                    *\n");
	out->printf("                           *     MINICHEM %s   *\n", MINICHEM_VERSION);
	out->printf("                           *                    *\n");
	out->printf("                           **********************\n");
	out->println();
	out->printf("Minichem is an open-source simple quantum chemistry code, \
written in educational purposes\n");
	out->printf("Build date: %s %s\n", MINICHEM_BUILD_DATE, MINICHEM_BUILD_TIME);
	#if defined __ICC
	out->printf("Compiler:   Intel C Compiler %d\n", __ICC);
	#elif defined __GNUC__
	out->printf("Compiler:   gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
	#else
	out->printf("Compiler:   undetected\n");
	#endif
	printf("Date:       %s", asctime(localtime(&t)));
	out->printf("Authors:\n");
	out->printf("  Alexander Oleynichenko        alexvoleynichenko@gmail.com\n");
	out->println();
}

int main(int argc, char **argv)
{
	// init
	startTimePoint = std::chrono::high_resolution_clock::now();
	mainlog = new Log("minichem.log");

	// run minichem's kernel with command-line parameters
	outputHeader();
	Kernel kernel(argc, argv);
	try {
		kernel.start();
	} catch (exception& e) {
		OutputStream::getStderr()->printf("Fatal error in minichem kernel: %s\n", e.what());
	}

	delete mainlog;

	return 0;
}






