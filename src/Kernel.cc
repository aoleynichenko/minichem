#include <chrono>
#include <stdexcept>

#include "Kernel.h"

using std::runtime_error;

namespace minichem {

Kernel::Kernel(int argc, char **argv)
{
	mainlog->log("Kernel::Kernel(int argc, char **argv) invoked");
}

int Kernel::start()
{
	mainlog->log("Starting kernel (Kernel::start())");
	return 0; // OK!
}

Kernel::~Kernel()
{

}

double Kernel::secondsFromStart()
{
	auto currTime = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::seconds>(currTime - startTimePoint).count();
	//cout << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << ":";
}

} // namespace minichem
