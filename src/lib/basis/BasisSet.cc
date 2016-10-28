#include <cctype>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "BasisSet.h"
#include "../chem/Elements.h"
#include "../chem/Molecule.h"
#include "../qscript/QS_Object.h"

namespace minichem {

using std::endl;
using std::fixed;
using std::invalid_argument;
using std::ostringstream;
using std::setprecision;
using std::set;
using std::setw;
using std::string;

BasisSet::BasisSet()
  : qscript::QS_Object()
{
  type = TYPE_BAS;
}

void BasisSet::addLBlock(std::string elemSym, LBlock block)
{
  // firstly, unify symbol of the chemical element
  elemSym[0] = toupper(elemSym[0]);
  if (elemSym.length() == 2)
    elemSym[1] = tolower(elemSym[1]);
  // add new block
  int Z = Elements::sym2charge(elemSym);
  set_[Z].push_back(block);
}

BasisSet BasisSet::filter(Molecule* mol) const
{
  BasisSet newbs;
  set<int> uniqElems = mol->uniqueElems();
  for (auto kv : set_)
    if (uniqElems.find(kv.first) != uniqElems.end())
      newbs.set_[kv.first] = kv.second;
  return newbs;
}

// I'm a great fan of NWChem... So, NWChem style!
string BasisSet::toString() const
{
  ostringstream str;  // stream in which we will output basis set data

  for (auto it = set_.begin(); it != set_.end(); it++) {
    int Z = it->first;
    BasisTemplate_t bt = it->second;
    Element elem = Elements::getElementByZ(Z);
    // print header
    string head = elem.sym + " (" + elem.name + ")";
    str << "  " << head << endl << "  ";
    for (int c = 0; c < head.length(); c++)
      str << "-";
    str << endl;
    // print functions
    int count = 1; // enumerate shells of element 'elem'
    for (auto block = bt.begin(); block != bt.end(); block++)
      for (auto c = block->contr_.begin(); c != block->contr_.end(); c++) {
        str << setw(4) << count << setw(4) << BasisSet::am2string(block->l_);
        for (auto prim = 0; prim < c->size(); prim++)
          str << fixed << setprecision(8) << setw(18+8*(prim!=0))
              << block->alpha_[prim] << setw(15) << (*c)[prim] << endl;
        count++;
      }
    str << endl;
  }

  return str.str();  // str is an ostringstream object!
}

string BasisSet::getTypeString() const
{
  return "basisset";
}

static string str_tolower(string s)
{
	string t = s;
	for (unsigned int i = 0; i < s.length(); i++)
		t[i] = tolower(s[i]);
	return t;
}

string BasisSet::am2string(int L)  // L -> {S, P, D, F...}
{
  if (L < 0 || L > 4)
    throw invalid_argument("wrong angular momentum integer: should be in range [0:4]");
  static string momlabels[] = {"S", "P", "D", "F", "G"};
  return momlabels[L];
}

int BasisSet::parseAngmom(std::string ams)  // {S, P, D, F... -> L}
{
  ams = str_tolower(ams);
	static string momlabels[] = {"s", "p", "d", "f", "g"};
	int maxmom = 4;
	for (int i = 0; i < maxmom; i++)
		if (momlabels[i] == ams)
			return i;
	return -1;
}

} //namespace minichem
