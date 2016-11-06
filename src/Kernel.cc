#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "Kernel.h"
#include "./lib/basis/BasisSet.h"
#include "./lib/chem/Elements.h"
#include "./lib/chem/Molecule.h"
#include "./lib/except/SyntaxError.h"
#include "./lib/io/Lexer.h"
#include "./lib/io/OutputStream.h"
#include "./lib/io/Token.h"

// QuantumScript engine essentials
#include "./lib/qscript/QS_Object.h"
#include "./lib/qscript/QS_RuntimeError.h"
#include "./lib/qscript/QS_Scope.h"

// Quantum-chemical subroutines
#include "./scf/scf.h"

using std::exception;
using std::ifstream;
using std::ostringstream;
using std::runtime_error;
using std::shared_ptr;
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

	// find basis set libraries using environmental variable MINICHEM_HOME
	char* homedir = getenv("MINICHEM_HOME");
	if (homedir) {
		libPath = string(homedir);
		mainlog->log("MINICHEM_HOME = " + libPath);
	}
	else
		mainlog->log("MINICHEM_HOME variable is void, you can set it to simplify \
search and inclusion of basis sets");

	// init QuantumScript engine
	currBasis = nullptr;
	currMolecule = nullptr;

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

	auto t0 = std::chrono::system_clock::now();

	// firstly, read settings from ~/.minichemrc
	// .minichemrc is an ordinary minichem script, user can write there,
	// for example, path to basis libraries, scf options (bad style!!!), etc.
	// TROUBLE WITH TILDA :(
	/*
	if (fileExists("~/.minichemrc")) {
		mainlog->log(".minichemrc file found, trying to read settings");
		ifstream rc("~/.minichemrc");
		lex.setInput(&rc);
		try {
			execScript();
		} catch (exception& e) {
			mainlog->log("[ERROR] While reading ~/.minichemrc file: %s",   e.what());
			out->printf ("[ERROR] While reading ~/.minichemrc file: %s\n", e.what());
		}
	}*/

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

	// print wall timing
	// this code may contain errors :)
	auto t1 = std::chrono::system_clock::now();
	double wallt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	int ms  = ((int) wallt) % 1000;
	int sec = ((int) wallt - ms)/1000;
	int min = (sec - sec % 60) / 60;
	int hrs = (min - min % 60) / 60;
	int day = (hrs - hrs % 24) / 24;
	out->printf("Wall time: %d days %d hours %d min %d sec %d ms\n", day, hrs%24, min%60, sec%60, ms);

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

OutputStream* Kernel::getOutput()
{
	return out;
}

inline void Kernel::hline()
{
	out->printf("---------------------------------------------------\
-----------------------------\n");   // 80 x '-'
}

/****************************************************************************/
/***                                PARSER                                 ***/
/*****************************************************************************/

void Kernel::execScript()
{
	Token t = lex.get();
	while (t.ttype != Token::TT_EOF) {
		if (t.ttype == Token::TT_KW_MOL)
			declMolecule();
		else if (t.ttype == Token::TT_KW_BASIS) {
			BasisSet* bs = new BasisSet();
			lex.putback(t);  // we will start basis set declaration from [basis] token
			declBasisSet(bs, &lex);
			scope_m.set("bs", bs);
		}
		else if (t.ttype == Token::TT_KW_TASK)
			runTask();
		else if (t.ttype == Token::TT_KW_TYPEOF)
			doTypeof();
		else if (t.ttype == Token::TT_KW_PRINT)
			doPrint();
		else if (t.ttype == Token::TT_KW_CURR)
			doCurr();
		else {
			mainlog->log("[ERROR] Unknown token in Kernel::execScript(): " + t.toString());
			throw SyntaxError("unknown token: " + t.toString());
		}
		t = lex.get();
	}
	mainlog->log("End of input file (EOF)");
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
	int bohrs = 1;  // bohrs by default
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
	Molecule* mol = new Molecule();
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
		mol->addAtom(Z, c*x, c*y, c*z);
		t = lex.get();
	}

	mol->setMult(mult);
	mol->setCharge(charge);

	if (t.ttype != '}') {
		mainlog->log("[ERROR] In Kernel::declMolecule(): '}' expected");
		throw SyntaxError("expected '}' in molecule declaration");
	}
	mol->check();  // check nelec, charge and multiplicity
	// all is OK!
	scope_m.set(molname, mol);
	mainlog->log("Molecule specification was succesfully read, charge = %d, \
mult = %d, nelec = %d", mol->getCharge(), mol->getMult(), mol->nelec());
}

// this function is designed to be very flexible and recursive in order to read
// basis sets from nested 'basis' instructions
void Kernel::declBasisSet(BasisSet* bs, Lexer* lexer)
{
	mainlog->log("Basis set declaration started");

	Token t = lexer->get(); // skip [KEYWORD|basis] token
	t = lexer->get();
	string setname = "unnamed";
	if (t.ttype == Token::TT_WORD || t.ttype == Token::TT_QUOTE) {
		setname = t.sval;
		t = lexer->get();
	}
	if (t.ttype != '{') {
		if (setname != "unnamed") {  // example: basis "cc-pvtz"
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
			// 3. search in minichem's basis library
			else if (fileExists(this->libPath + "/lib/basissets/" + setname)) {
				mainlog->log("Search basis set '%s' in minichem's library", setname.c_str());
				string path = this->libPath + "/lib/basissets/" + setname;
				ifstream subbas(path);
				Lexer sublex;
				sublex.setInput(&subbas);
				mainlog->log("File '%s' exists (in library directory), invoke recursive", setname.c_str());
				declBasisSet(bs, &sublex);
			}
			else
				throw SyntaxError("basis set '" + setname + "' not found");
			lex.putback(t);  // because se have tested input for {
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
				throw SyntaxError("wrong angular momentum: " + am.sval +
					" (expected S, P, D, F or G, found '" + am.sval + "')");
			// read block!
			lexer->setEolEnabled(true);
			BasisSet::LBlock block;
			block.l_ = L;
			block.cart_ = cartesian;
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
				vector<double> coeffs;  // one HORIZONTAL line with coefficients
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
 				if (line_n == 1) {
					ncontr = coeffs.size();
					for (int i = 0; i < ncontr; i++)
						block.contr_.push_back(vector<double>()); // create NCONTR contracted functions
				}
				else  // line_n > 1
					if (ncontr != coeffs.size()) {
						throw SyntaxError("in basis set declaration: expected rectangular \
matrix of contraction coefficients");}
				// all is OK
				for (int i = 0; i < ncontr; i++)
					block.contr_[i].push_back(coeffs[i]);  // horizontal vector -> vectical vectors
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
			bs->addLBlock(sym, block);
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

void Kernel::doTypeof()
{
	vector<string> names;

	lex.setEolEnabled(true);
	Token t(Token::TT_EOF);
	do {
		t = lex.get();
		if (t.ttype != Token::TT_WORD)
			throw SyntaxError("in typeof operator: expected word, but found " + t.toString());
		names.push_back(t.sval);
		t = lex.get();
	} while (t.ttype == ',');
	if (t.ttype != Token::TT_EOL && t.ttype != Token::TT_EOF)
		throw SyntaxError("expected end-of-line after 'typeof' statement");
	lex.setEolEnabled(false);

	for (auto name : names)
		out->printf("%s\n", scope_m.get(name)->getTypeString().c_str());
}

void Kernel::doPrint()
{
	vector<string> names;

	lex.setEolEnabled(true);
	Token t(Token::TT_EOF);
	do {
		t = lex.get();
		if (t.ttype != Token::TT_WORD)
			throw SyntaxError("in 'print' operator: expected word, but found " + t.toString());
		names.push_back(t.sval);
		t = lex.get();
	} while (t.ttype == ',');
	if (t.ttype != Token::TT_EOL && t.ttype != Token::TT_EOF)
		throw SyntaxError("expected end-of-line after 'print' statement");
	lex.setEolEnabled(false);

	for (auto name : names)
		out->printf("%s\n", scope_m.get(name)->toString().c_str());
}

void Kernel::doCurr()
{
	using namespace qscript;
	vector<string> names;

	lex.setEolEnabled(true);
	Token t(Token::TT_EOF);
	do {
		t = lex.get();
		if (t.ttype != Token::TT_WORD)
			throw SyntaxError("in 'curr' operator: expected word, but found " + t.toString());
		names.push_back(t.sval);
		t = lex.get();
	} while (t.ttype == ',');
	if (t.ttype != Token::TT_EOL && t.ttype != Token::TT_EOF)
		throw SyntaxError("expected end-of-line after 'print' statement");
	lex.setEolEnabled(false);

	for (auto name : names) {
		QS_Object* obj = scope_m.get(name);
		if (obj->type == QS_Object::TYPE_BAS)
			currBasis = (BasisSet*) obj;
		else if (obj->type == QS_Object::TYPE_MOL)
			currMolecule = (Molecule*) obj;
	}
}

void Kernel::runTask()
{
	Token t = lex.get();
	if (t.ttype != Token::TT_WORD)
		throw SyntaxError("in task directive: expected one of [scf], but found " + t.toString());
	if (t.sval == "scf") {  // run Hartree-Fock calculation
		if (!currMolecule || !currBasis)
			throw runtime_error("please, specify current molecule and basis set");
		// auto selection of proper HF method
		// Note!!! You may have singlet with na != nb (of course, it is not ground state).
		// Setting number of electrons for orbitals in each irrep is much better!
		// (to be implemented in future)
		auto na = currMolecule->nalpha();
		auto nb = currMolecule->nbeta();
		if (na == nb)
			rhf(this, currBasis, currMolecule);
		else
			uhf(this, currBasis, currMolecule);
	}
	else
		throw SyntaxError("in task directive: " + t.sval + " method is not yet implemented'");
}

// helper functions
std::string getOsName()
{
	string bits = "";
	#ifdef __i386__
	bits = " x86";
	#endif
	#ifdef __x86_64__
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

static string str_tolower(string s)
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
	try {
		return Elements::sym2charge(s);
	}
	catch (exception&) {
		return 0;  // symbol not found
	}
}

int parseAngularMomentum(string ams)
{
	ams = str_tolower(ams);
	string momlabels[] = {"s", "p", "d", "f", "g"};
	int maxmom = 5;
	for (int i = 0; i < maxmom; i++)
		if (momlabels[i] == ams)
			return i;
	return -1;
}

} // namespace minichem
