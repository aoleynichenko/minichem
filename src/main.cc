/* Minichem - open-source simple quantum chemistry code,
 * written in eduactional purposes.
 * 
 * Alexander Oleynichenko, Moscow State University, 2016
 * alexvoleynichenko@gmail.com
 *
 *  main.cc
 *  Kernel of minichem application.
 *  Invokes all another modules: input, hf, mp2, cc, geom, etc.
 */

#include "OutputStream.h"

using namespace minichem;

int main(int argc, char **argv)
{
	OutputStream::getStdout()->printf("Hello, Quantum World!\n");
	OutputStream::getStderr()->printf("Errors :(\n");
	return 0;
}







