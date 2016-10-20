#ifndef _MOLECULE_H_INCLUDED
#define _MOLECULE_H_INCLUDED

#include <string>
#include <vector>

namespace minichem {

class Molecule {
public:
	struct Atom {
		double x, y, z;
		int charge;
		Atom(int, double, double, double);
	};

	Molecule();
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
	std::string toString();
	void check() const;
private:
	int mult;
	int charge;
	std::vector<Atom> atoms;

	struct Element {
		std::string sym;
		int Z;
		double mass;
		Element(int z, std::string s, double m);
	};
	struct PeriodicTable {
		std::vector<Element> elements;
		PeriodicTable();
		void addElement(Element el);
		Element& getElementByZ(int z);
		Element& getElementBySym(std::string sym);
	};
	static PeriodicTable perTable;
};

typedef Molecule::Atom Atom;

} // namespace minichem

#endif // _MOLECULE_H_INCLUDED