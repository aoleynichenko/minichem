#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "Kernel.h"
#include "./lib/chem/Molecule.h"
#include "./lib/except/SyntaxError.h"
#include "./lib/io/Lexer.h"
#include "./lib/io/OutputStream.h"
#include "./lib/io/Token.h"

using std::ifstream;
using std::ostringstream;
using std::runtime_error;
using std::string;
using std::vector;

namespace minichem {

std::string getOsName();
bool fileExists(const std::string& name);
int  isElementSymbol(string s);
bool isint(double d);

Kernel::Kernel(int argc, char **argv)
{
	mainlog->log("Kernel::Kernel(int argc, char **argv) invoked, argc == %d", argc);
	char hostname[128];
	gethostname(hostname, 128);
	mainlog->log("Running on %s (%s)", hostname, getOsName().c_str());
	ostringstream params;
	for (int i = 0; i < argc; i++)
		params << argv[i] << " ";
	mainlog->log("Kernel arguments: \"%s\"", params.str().c_str());

	// attach "standard" output
	out = OutputStream::getStdout();

	// which input files should we process?
	// inputFiles_m is implicitly initialized
	for (int i = 1; i < argc; i++)
		inputFiles_m.push_back(string(argv[i]));
	params.str("");
	params.clear();
	for (auto inp = inputFiles_m.begin(); inp != inputFiles_m.end(); inp++) {
		params << *inp << " ";
		if (!fileExists(*inp)) {
			mainlog->log("[WARNING] File \"%s\" doesn't exists!", inp->c_str());
			out->printf("[WARNING] File \"%s\" doesn't exists!\n", inp->c_str());
		}
	}
	mainlog->log("Input files: %d ( %s)", inputFiles_m.size(), params.str().c_str());
}

int Kernel::start()
{
	mainlog->log("Starting kernel (Kernel::start())");

	for (auto inp = inputFiles_m.begin(); inp != inputFiles_m.end(); inp++) {
		hline();
		if (!fileExists(*inp)) {
			mainlog->log("[ERROR] File \"%s\" doesn't exists, will be ignored", inp->c_str());
			out->printf("[ERROR] File \"%s\" doesn't exists, will be ignored\n", inp->c_str());
			continue; // no fatal error
		}
		out->printf("Input file: %s\n", inp->c_str());
		mainlog->log("New input file: \"%s\"", inp->c_str());
		// open and execute input file line-by-line
		ifstream input(*inp);
		lex.setInput(&input);
		try {
			execScript();
		} catch (SyntaxError& se) {
			mainlog->log("[ERROR] Syntax error in file %s: %s", inp->c_str(), se.what());
			out->printf ("[ERROR] Syntax error in file %s: %s\n", inp->c_str(), se.what());
		} catch (runtime_error& re) {
			mainlog->log("[ERROR] Runtime error in file %s: %s", inp->c_str(), re.what());
			out->printf ("[ERROR] Runtime error in file %s: %s\n", inp->c_str(), re.what());
		}
	}


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

inline void Kernel::hline()
{
	out->printf("---------------------------------------------------\
-----------------------------\n");   // 80 x '-'
}

/*****************************************************************************/
/***                                PARSER                                 ***/
/*****************************************************************************/

void Kernel::execScript()
{
	Token t = lex.get();
	if (t.ttype == Token::TT_KW_MOL)
		declMolecule();
	else {
		mainlog->log("[ERROR] Unknown token in Kernel::execScript(): " + t.toString());
		throw SyntaxError("unknown token: " + t.toString());
	}
}

void Kernel::declMolecule()
{
	mainlog->log("Molecule declaration started");

	Token t = lex.get();
	string molname = "default";
	if (t.ttype == Token::TT_WORD) {
		molname = t.sval;
		t = lex.get();
	}
	if (t.ttype != '{') {
		mainlog->log("[ERROR] In Kernel::declMolecule(): '{' expected");
		throw SyntaxError("expected '{' in molecule declaration");
	}

	// get additional parameters
	int mult = 1;
	int charge = 0;
	int bohrs = 1;  // angstroms by default
	t = lex.get();
	while (t.ttype == Token::TT_WORD) {
		if (isElementSymbol(t.sval))
			break;
		if (t.sval == "mult") {
			mult = lex.getint();
			if (mult <= 0)
				throw SyntaxError("multiplicity should be strictly positive");
		}
		else if (t.sval == "charge")
			charge = lex.getint();
		else if (t.sval == "units") {
			t = lex.get();
			if (t.ttype != Token::TT_WORD || (t.sval != "atomic" && t.sval != "angstroms"))
				throw SyntaxError("after 'units' keyword: 'atomic' or 'angstroms', found: "
					+ t.toString());
			if (t.sval == "angstroms")
				bohrs = 0;
		}
		else
			throw SyntaxError("unknown keyword in molecule declaration: \'" + t.sval + "\'");
		t = lex.get();
	}

	// read xyz and create new Molecule
	// now we have Symbol or integer Z
	Molecule mol;
	while (t.ttype == Token::TT_WORD || t.ttype == Token::TT_NUMBER) {
		int Z = 0;
		if (t.ttype == Token::TT_WORD) {  // sym  x  y  z
			if (!(Z = isElementSymbol(t.sval)))
				throw SyntaxError("is not a symbol of chemical element: \'" + t.sval + "\'");
		}
		else if (t.ttype == Token::TT_NUMBER) {  // Z  x  y  z
			if (!isint(t.dval))
				throw SyntaxError("charge of atom should be integer");
			Z = (int) t.dval;
			if (Z <= 0 || Z > 18)
				throw SyntaxError("wrong atomic charge, should be 0 < Z < 18 (H ... Ar)");
		}
		double x = lex.getdouble();
		double y = lex.getdouble();
		double z = lex.getdouble();
		double c = bohrs ? 1.0 : 1.889725989;  // convert A to a.u.
		mol.addAtom(Z, c*x, c*y, c*z);
		t = lex.get();
	}

	mol.setMult(mult);
	mol.setCharge(charge);

	if (t.ttype != '}') {
		mainlog->log("[ERROR] In Kernel::declMolecule(): '}' expected");
		throw SyntaxError("expected '}' in molecule declaration");
	}
	mol.check();  // check nelec, charge and multiplicity
	// all is OK!
	mainlog->log("Molecule specification was succesfully read, charge = %d, \
mult = %d, nelec = %d", mol.getCharge(), mol.getMult(), mol.nelec());
}


// helper functions
std::string getOsName()
{
	string bits = "";
	#ifdef __i386__
	bits = " x86";
	#endif
	#ifdef __x8_64__
	bits = " x86_64";
	#endif

    #ifdef _WIN32
    return "Windows 32-bit";
    #elif _WIN64
    return "Windows 64-bit";
    #elif __linux__
    return "Linux" + bits;
    #elif __APPLE__ || __MACH__
    return "Mac OSX" + bits;
    #elif __FreeBSD__
    return "FreeBSD" + bits;
    #elif __unix || __unix__
    return "Unix" + bits;
    #else
    return "Other";
    #endif
} 

inline bool fileExists(const std::string& name) {
    if (FILE* file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}

string str_tolower(string s)
{
	string t = s;
	for (unsigned int i = 0; i < s.length(); i++)
		t[i] = tolower(s[i]);
	return t;
}

bool isint(double d)
{
  if (fabs(d) != fabs((int) d))
    return false;
  return true;
}

int isElementSymbol(string s)
{
	static const int lensym = 18;
	static string symbols[] = {
		"h",  "he",
		"li", "be", "b",  "c",  "n",  "o",  "f",  "ne",
		"na", "mg", "al", "si", "p",  "s",  "cl", "ar"
	};
	s = str_tolower(s);
	for (size_t i = 0; i < lensym; i++)
		if (s == symbols[i])
			return i + 1;
	return 0;  // symbol not found
}

} // namespace minichem
