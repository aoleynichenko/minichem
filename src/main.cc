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

#include "minichem.h"
#include "OutputStream.h"

using namespace minichem;

void outputHeader()
{
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
	out->printf("Built: %s %s\n", MINICHEM_BUILD_DATE, MINICHEM_BUILD_TIME);
	out->printf("Authors:\n");
	out->printf("  Alexander Oleynichenko        alexvoleynichenko@gmail.com\n");
	out->println();
}

int main(int argc, char **argv)
{
	outputHeader();
	return 0;
}







