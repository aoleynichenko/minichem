#include <cctype>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Elements.h"
#include "Molecule.h"
#include "../io/OutputStream.h"
#include "../qscript/QS_Object.h"

namespace minichem {

using std::fixed;
using std::invalid_argument;
using std::ostringstream;
using std::runtime_error;
using std::setprecision;
using std::setw;
using std::string;
using std::vector;

//Molecule::PeriodicTable Molecule::perTable = Molecule::PeriodicTable();

Molecule::Atom::Atom(int Z, double _x, double _y, double _z)
	:	x(_x), y(_y), z(_z), charge(Z)
{
}

Molecule::Molecule()
	:	qscript::QS_Object(), mult(1), charge(0), atoms(vector<Atom>())
{
	type = TYPE_MOL;
}

Molecule::~Molecule()
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

void Molecule::addAtom(std::string symbol, double x, double y, double z)
{
	int Z = Elements::getElementBySym(symbol).Z;  // atomic charge
	addAtom(Z, x, y, z);
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
		Element e = Elements::getElementByZ(a->charge);
		mass += e.mass;
	}
	return mass;
}

// to be implemented!
//int Molecule::nalpha() const;
//int Molecule::nbeta()  const;
int Molecule::nelec() const
{
	int ne = 0;
	for (auto a = atoms.begin(); a != atoms.end(); a++)
		ne += a->charge;
	return ne - this->charge;
}

double Molecule::enuc() const
{
	double en = 0.0;
	for (size_t i = 0; i < atoms.size(); i++)
		for (size_t j = i + 1; j < atoms.size(); j++) {
			double dx = atoms[i].x - atoms[j].x;
			double dy = atoms[i].y - atoms[j].y;
			double dz = atoms[i].z - atoms[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			en += atoms[i].charge * atoms[j].charge / r;
		}
	return en;
}

string Molecule::toString() const
{
	ostringstream out;
	out << "Charge = " << charge << "\nMult = " << mult << "\n";
	out << fixed << setprecision(8);
	for (auto atom : atoms)
		out << "  " << std::left << setw(3) << Elements::charge2sym(atom.charge) << std::right
			<< setw(15) << atom.x << setw(15) << atom.y << setw(15) << atom.z << "\n";
	return out.str();
}

string Molecule::getTypeString() const
{
  return "molecule";
}

void Molecule::check() const
{
	int ne = nelec();
	int nonpaired = mult-1;
	if ((ne - nonpaired) % 2 != 0) {
		ostringstream errmsg;
		errmsg << "illegal charge/multiplicity: nelec = "
			<< ne << ", non-paired = " << nonpaired;
		throw runtime_error(errmsg.str());
	}
}

std::set<int> Molecule::uniqueElems() const
{
	std::set<int> uniq;
	for (auto atom : atoms)
		uniq.insert(atom.charge);
	return uniq;
}

} //namespace minichem
