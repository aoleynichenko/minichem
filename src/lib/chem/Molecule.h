#ifndef _MOLECULE_H_INCLUDED
#define _MOLECULE_H_INCLUDED

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "../qscript/QS_Object.h"

namespace minichem {

class Molecule : public qscript::QS_Object {
public:
	struct Atom {
		double x, y, z;
		int charge;
		Atom(int, double, double, double);
	};

	Molecule();
	~Molecule();
	std::vector<Atom>* getAtoms();
	void addAtom(int Z, double x, double y, double z);
	void addAtom(std::string symbol, double x, double y, double z);
	void addAtom(Atom a);
	void setMult(int mult);
	int  getMult() const;
	void setCharge(int charge);
	int  getCharge() const;
	int  natoms() const;
	double mass() const;
	int  nalpha() const;
	int  nbeta()  const;
	int  nelec()  const;
	double enuc() const;
	std::string toString() const;
	std::string getTypeString() const;
	void check() const;
	std::set<int> uniqueElems() const;
private:
	int mult;
	int charge;
	std::vector<Atom> atoms;
};

typedef Molecule::Atom Atom;
typedef std::shared_ptr<Molecule> Mol_ptr;

} // namespace minichem

#endif // _MOLECULE_H_INCLUDED
