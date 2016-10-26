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
#include "./lib/basis/BasisSet.h"
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
int parseAngularMomentum(string ams);
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
	else if (t.ttype == Token::TT_KW_BASIS) {
		BasisSet bs;
		lex.putback(t);  // we will start basis set declaration from [basis] token
		declBasisSet(&bs, &lex);
	}
	else {
		mainlog->log("[ERROR] Unknown token in Kernel::execScript(): " + t.toString());
		throw SyntaxError("unknown token: " + t.toString());
	}
}

void Kernel::declMolecule()
{
	mainlog->log("Molecule declaration started");

	Token t = lex.get();
	string molname = "unnamed";
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

// this function is designed to be very flexible and recursive in order to read
// basis sets from nested 'basis' instructions
void Kernel::declBasisSet(BasisSet* bs, Lexer* lexer)
{
	mainlog->log("Basis set declaration started");

	Token t = lexer->get(); // skip [KEYWORD|basis] token
	t = lexer->get();
	string setname = "unnamed";
	if (t.ttype == Token::TT_WORD) {
		setname = t.sval;
		t = lexer->get();
	}
	if (t.ttype != '{') {
		if (setname != "unnamed") {  // basis cc-pvtz
			mainlog->log("In Kernel::declBasisSet(): basis set '%s' declared in short notation", setname.c_str());
			// load basis from file/library or make it current if it is already loaded
			// 1. search in our variables
			// 2. search in current directory
			if (fileExists(setname)) {
				ifstream subbas(setname);
				Lexer sublex;
				sublex.setInput(&subbas);
				mainlog->log("File '%s' exists, invoke recursive", setname.c_str());
				declBasisSet(bs, &sublex);
			}
			else {
				// 3. search in minichem's basis library
				mainlog->log("Search basis set '%s' in minichem's library", setname.c_str());
				throw SyntaxError("basis set '" + setname + "' not found");
			}
			return;
		}
		else {
			mainlog->log("[ERROR] In Kernel::declBasisSet(): '{' expected, but found " + t.toString());
			throw SyntaxError("expected '{' in basis set declaration, but found " + t.toString());
		}
	}

	// loop over keywords and element symbols
	// * -> invoke library
	// Sym -> invoke library OR read L-block (block with definite angular momentum)
	// Keyword: one of:
	//  cartesian
	//  spherical
	bool cartesian = false; // by default, spherical basis sets
	t = lexer->get();
	while (t.ttype == Token::TT_WORD) {
		if (t.sval == "cartesian")
			cartesian = true;
		else if (t.sval == "spherical")
			cartesian = false;
		else if (isElementSymbol(t.sval)) {  // read L-block
			string sym = t.sval;  // save element symbol
			Token am = lexer->get();
			if (am.ttype != Token::TT_WORD)
				throw SyntaxError("wrong angular momentum should be string, but found " + am.toString());
			int L = parseAngularMomentum(am.sval);
			if (L < 0)
				throw SyntaxError("wrong angular momentum: " + am.sval + " (expected S, P, D, F ot G)");
			// read block!
			lexer->setEolEnabled(true);
			BasisSet::LBlock block;
			block.l_ = L;
			int line_n = 1;
			int ncontr = 0;  // number of contraction coeffs in line, should be the same for all lines
			t = lexer->get();
			t = lexer->match(t, Token::TT_EOL);  // End-Of-Line after angular momentum symbol!
			// block ending conditions:
			//   '{'
			//   Sym
			// alpha may be a string (if variable)
			while (t.ttype == Token::TT_NUMBER ||
				(t.ttype == Token::TT_WORD && !isElementSymbol(t.sval))) {
				// token == alpha
				// no variables in minichem 0.1
				// parse ONE string: <alpha> <c1> <c2> ... <cN>
				if (t.ttype == Token::TT_NUMBER) {
					block.alpha_.push_back(t.dval);
				}
				else {
					throw SyntaxError("in basis set declaration: variables are not allowed\
yet, alpha should be a number, but found " + t.toString());
				}
				vector<double> coeffs;
				t = lexer->get();
				while (t.ttype != Token::TT_EOL && t.ttype != Token::TT_EOF) {
					if (t.ttype != Token::TT_NUMBER)
						throw SyntaxError("in basis set declaration: variables are not allowed \
yet, contraction coefficient should be a number, but found " + t.toString());
					coeffs.push_back(t.dval);
					t = lexer->get();
				}
				while (t.ttype == Token::TT_EOL) // skip white space between lines in block
					t = lexer->get();
				if (t.ttype == Token::TT_EOF)
					throw SyntaxError("in basis set declaration: unexpected end of file");
				lexer->putback(t);  // return next token (!= EOF!), maybe it is alpha

				// verify line of contraction coefficient
				if (line_n == 1 && coeffs.size() == 0)
					throw SyntaxError("in basis set declaration: the number of contraction\
 coefficients should be non-zero");
 				if (line_n == 1)
					ncontr = coeffs.size();
				else  // line_n > 1
					if (ncontr != coeffs.size()) {
						throw SyntaxError("in basis set declaration: expected rectangular \
matrix of contraction coefficients");}
				// all is OK
				block.contr_.push_back(coeffs);
				line_n++;
				// get next alpha, right curly bracket or Sym
				t = lexer->get();
			}
			lexer->putback(t);
			lexer->setEolEnabled(false);
			// now we have valid block ("template") of contracted GTOs
			// write info to log for debugging and add this block to BasisSet object
			ostringstream alphas;  // too complicated:) C++, where is join()?
			for (size_t i = 0; i < block.alpha_.size(); i++) {
				alphas << block.alpha_[i];
				if (i != block.alpha_.size()-1)
					alphas << ",";
			}
			mainlog->log("Succesfully read new L-block for element %s: {L=%d, alpha=[%s], \
Ncontracted=%d}", sym.c_str(), block.l_, alphas.str().c_str(), ncontr);
		}
		else {
			mainlog->log("[ERROR] In Kernel::declBasisSet(): is not an element symbol: '%s'", t.sval.c_str());
			throw SyntaxError("in basis set declaration: is not an element symbol: " + t.sval);
		}
		t = lexer->get();
	}
	if (t.ttype != '}') {
		mainlog->log("[ERROR] In Kernel::declBasisSet(): '}' expected, but found " + t.toString());
		throw SyntaxError("expected '}' in basis set declaration");
	}
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

int parseAngularMomentum(string ams)
{
	ams = str_tolower(ams);
	string momlabels[] = {"s", "p", "d", "f", "g"};
	int maxmom = 4;
	for (int i = 0; i < maxmom; i++)
		if (momlabels[i] == ams)
			return i;
	return -1;
}

} // namespace minichem
