/* Kernel.h
 * This header contains definition of class Kernel - main class
 * of the Minichem program, which contains all runtime environment - 
 * variables, settings, results of previous calculations, memory pool.
 *
 * Note: maybe, the kernel can be stored as a file?
 *
 * Author: A. Oleynichenko, 2016.
 * Mailto: alexvoleynichenko@gmail.com
 */

#ifndef _KERNEL_H_INCLUDED
#define _KERNEL_H_INCLUDED

#include <chrono>

#include "./lib/io/Log.h"

namespace minichem {

class Kernel {
public:
	Kernel(int argc, char **argv);
	int start();

	~Kernel();

	static double secondsFromStart();
};

// Main log. It is separated from kernel variables to be stable wrt kernel crash.
// For initialization, see main.cc.
extern Log* mainlog;
extern std::chrono::high_resolution_clock::time_point startTimePoint;

} // namespace minichem

#endif // _KERNEL_H_INCLUDED