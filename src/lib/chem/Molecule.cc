#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Molecule.h"

namespace minichem {

using std::invalid_argument;
using std::ostringstream;
using std::string;
using std::vector;

Molecule::PeriodicTable Molecule::perTable = Molecule::PeriodicTable();

Molecule::Atom::Atom(int Z, double _x, double _y, double _z)
	:	x(_x), y(_y), z(_z), charge(Z)
{
}

Molecule::Molecule()
	:	mult(1), charge(0), atoms(vector<Atom>())
{
}

void Molecule::addAtom(Atom a)
{
	atoms.push_back(a);
}

void Molecule::addAtom(int Z, double x, double y, double z)
{
	atoms.push_back(Atom(Z, x, y, z));
}

vector<Atom>* Molecule::getAtoms()
{
	return &atoms;
}

void Molecule::setMult(int _mult)
{
	mult = _mult;
}

int Molecule::getMult() const
{
	return mult;
}

void Molecule::setCharge(int _charge)
{
	charge = _charge;
}

int Molecule::getCharge() const
{
	return charge;
}

int Molecule::natoms() const
{
	return atoms.size();
}

double Molecule::mass() const
{
	double mass = 0.0;
	for (auto a = atoms.begin(); a != atoms.end(); a++) {
		Element e = perTable.getElementByZ(a->charge);
		mass += e.mass;
	}
	return mass;
}

// to be implemented!
//int Molecule::nalpha() const;
//int Molecule::nbeta()  const;
//int Molecule::nelec()  const;

string Molecule::toString()
{
	return "";
}

void Molecule::check() const
{

}

// private class Molecule::Element
Molecule::Element::Element(int z, string s, double m)
	:	 sym(s), Z(z), mass(m)
{
}


// private class Molecule::PeriodicTable
Molecule::PeriodicTable::PeriodicTable()
	:	elements(vector<Molecule::Element>())
{
	addElement(Element(1,  "H",  1));
	addElement(Element(2,  "He", 4));
	addElement(Element(3,  "Li", 7));
	addElement(Element(4,  "Be", 9));
	addElement(Element(5,  "B",  11));
	addElement(Element(6,  "C",  12));
	addElement(Element(7,  "N",  14));
	addElement(Element(8,  "O",  16));
	addElement(Element(9,  "F",  19));
	addElement(Element(10, "Ne", 20));
}

void Molecule::PeriodicTable::addElement(Element el)
{
	elements.push_back(el);
}

Molecule::Element& Molecule::PeriodicTable::getElementByZ(int z)
{
	for (auto p = elements.begin(); p != elements.end(); p++)
		if (p->Z == z)
			return *p;

	std::ostringstream errmsg;
	errmsg << "wrong atomic charge: " << z << " (not found in Mendeleev's table)";
	throw invalid_argument(errmsg.str());
}

string stolower(string s)
{
	string t = s;
	for (unsigned int i = 0; i < s.length(); i++)
		t[i] = tolower(s[i]);
	return t;
}

Molecule::Element& Molecule::PeriodicTable::getElementBySym(string sym)
{
	string s = stolower(sym);
	for (auto p = elements.begin(); p != elements.end(); p++)
		if (stolower(p->sym) == s)
			return *p;

	std::ostringstream errmsg;
	errmsg << "wrong symbol of element: " << sym << " (not found in Mendeleev's table)";
	throw invalid_argument(errmsg.str());
}

} //namespace minichem
