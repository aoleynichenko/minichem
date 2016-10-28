#include "Elements.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace minichem {

using std::invalid_argument;
using std::ostringstream;
using std::string;
using std::vector;

// idea from libint2/chemistry/elements.h
static std::vector<Element> info {
  {1,   1.0079 ,  "H" , "Hydrogen"     },
  {2,   4.0026 ,  "He", "Helium"       },
  {3,   6.941  ,  "Li", "Lithium"      },
  {4,   9.0122 ,  "Be", "Beryllium"    },
  {5,   10.811 ,  "B" , "Boron"        },
  {6,   12.0107,  "C" , "Carbon"       },
  {7,   14.0067,  "N" , "Nitrogen"     },
  {8,   15.9994,  "O" , "Oxygen"       },
  {9,   18.9984,  "F" , "Fluorine"     },
  {10,  20.1797,  "Ne", "Neon"         },
  {11,  22.9897,  "Na", "Sodium"       },
  {12,  24.305 ,  "Mg", "Magnesium"    },
  {13,  26.9815,  "Al", "Aluminum"     },
  {14,  28.0855,  "Si", "Silicon"      },
  {15,  30.9738,  "P" , "Phosphorus"   },
  {16,  32.065 ,  "S" , "Sulfur"       },
  {17,  35.453 ,  "Cl", "Chlorine"     },
  {18,  39.0983,  "Ar", "Argon"        },
  {19,  39.948 ,  "K" , "Potassium"    },
  {20,  40.078 ,  "Ca", "Calcium"      },
  {21,  44.9559,  "Sc", "Scandium"     },
  {22,  47.867 ,  "Ti", "Titanium"     },
  {23,  50.9415,  "V" , "Vanadium"     },
  {24,  51.9961,  "Cr", "Chromium"     },
  {25,  54.938 ,  "Mn", "Manganese"    },
  {26,  55.845 ,  "Fe", "Iron"         },
  {27,  58.6934,  "Co", "Cobalt"       },
  {28,  58.9332,  "Ni", "Nickel"       },
  {29,  63.546  , "Cu", "Copper"       },
  {30,  65.39   , "Zn", "Zinc"         },
  {31,  69.723  , "Ga", "Gallium"      },
  {32,  72.64   , "Ge", "Germanium"    },
  {33,  74.9216 , "As", "Arsenic"      },
  {34,  78.96   , "Se", "Selenium"     },
  {35,  79.904  , "Br", "Bromine"      },
  {36,  83.8    , "Kr", "Krypton"      },
  {37,  85.4678 , "Rb", "Rubidium"     },
  {38,  87.62   , "Sr", "Strontium"    },
  {39,  88.9059 , "Y" , "Yttrium"      },
  {40,  91.224  , "Zr", "Zirconium"    },
  {41,  92.9064 , "Nb", "Niobium"      },
  {42,  95.94   , "Mo", "Molybdenum"   },
  {43,  98      , "Tc", "Technetium"   },
  {44,  101.07  , "Ru", "Ruthenium"    },
  {45,  102.9055, "Rh", "Rhodium"      },
  {46,  106.42  , "Pd", "Palladium"    },
  {47,  107.8682, "Ag", "Silver"       },
  {48,  112.411 , "Cd", "Cadminium"    },
  {49,  114.818 , "In", "Indium"       },
  {50,  118.71  , "Sn", "Tin"          },
  {51,  121.76  , "Sb", "Antimony"     },
  {52,  126.9045, "Te", "Tellurium"    },
  {53,  127.6   , "I" , "Iodine"       },
  {54,  131.293 , "Xe", "Xenon"        },
  {55,  132.9055, "Cs", "Cesium"       },
  {56,  137.327 , "Ba", "Barium"       },
  {57,  138.9055, "La", "Lanthanium"   },
  {58,  140.116 , "Ce", "Cerium"       },
  {59,  140.9077, "Pr", "Praseodymium" },
  {60,  144.24  , "Nd", "Neodymium"    },
  {61,  145     , "Pm", "Promethium"   },
  {62,  150.36  , "Sm", "Samarium"     },
  {63,  151.964 , "Eu", "Europium"     },
  {64,  157.25  , "Gd", "Gadolinium"   },
  {65,  158.9253, "Tb", "Terbium"      },
  {66,  162.5   , "Dy", "Dysprosium"   },
  {67,  164.9303, "Ho", "Holmium"      },
  {68,  167.259 , "Er", "Erbium"       },
  {69,  168.9342, "Tm", "Thulium"      },
  {70,  173.04  , "Yb", "Ytterbium"    },
  {71,  174.967 , "Lu", "Lutetium"     },
  {72,  178.49  , "Hf", "Hafnium"      },
  {73,  180.9479, "Ta", "Tantalum"     },
  {74,  183.84  , "W" , "Tungsten"     },
  {75,  186.207 , "Re", "Rhenium"      },
  {76,  190.23  , "Os", "Osmium"       },
  {77,  192.217 , "Ir", "Iridium"      },
  {78,  195.078 , "Pt", "Platinum"     },
  {79,  196.9665, "Au", "Gold"         },
  {80,  200.59,   "Hg", "Mercury"      },
  {81,  204.3833, "Tl", "Thallium"     },
  {82,  207.2,    "Pb", "Lead"         },
  {83,  208.9804, "Bi", "Bismuth"      },
  {84,  209,      "Po", "Polonium"     },
  {85,  210,      "At", "Astatine"     },
  {86,  222,      "Rn", "Radon"        },
  {87,  223,      "Fr", "Francium"     },
  {88,  226,      "Ra", "Radium"       },
  {89,  227,      "Ac", "Actinium"     },
  {90,  232.0381, "Th", "Thorium"      },
  {91,  231.0359, "Pa", "Protactinium" },
  {92,  238.0289, "U" , "Uranium"      },
  {93,  237,      "Np", "Neptunium"    },
  {94,  244,      "Pu", "Plutonium"    },
  {95,  243,      "Am", "Americium"    },
  {96,  247,      "Cm", "Curium"       },
  {97,  247,      "Bk", "Berkelium"    },
  {98,  251,      "Cf", "Californium"  },
  {99,  252,      "Es", "Einsteinum"   },
  {100, 257,      "Fm", "Fermium"      },
  {101, 258,      "Md", "Mendelevium"  },
  {102, 259,      "No", "Nobelium"     },
  {103, 262,      "Lr", "Lawrencium"   },
  {104, 261,      "Rf", "Rutherfordium"},
  {105, 262,      "Db", "Dubnium"      },
  {106, 266,      "Sg", "Seaborgium"   },
  {107, 264,      "Bh", "Bohrium"      },
  {108, 277,      "Hs", "Hassium"      },
  {109, 268,      "Mt", "Meitnerium"   },
  {110, 0,        "Ds", "Darmstadtium" },
  {111, 272,      "Rg", "Roentgenium"  },
  {112, 0,        "Cn", "Copernicium"  },
  {113, 0,        "Nh", "Nihonium"     },
  {114, 0,        "Fl", "Flerovium"    },
  {115, 0,        "Mc", "Moscovium"    },
  {116, 0,        "Lv", "Livermorium"  },
  {117, 0,        "Ts", "Tennessine"   },
  {118, 0,        "Og", "Oganesson"    }
};

Element& Elements::getElementByZ(int z)
{
	for (auto p = info.begin(); p != info.end(); p++)
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

Element& Elements::getElementBySym(string sym)
{
	string s = stolower(sym);
	for (auto p = info.begin(); p != info.end(); p++)
		if (stolower(p->sym) == s)
			return *p;

	std::ostringstream errmsg;
	errmsg << "wrong symbol of element: " << sym << " (not found in Mendeleev's table)";
	throw invalid_argument(errmsg.str());
}

int Elements::sym2charge(string sym)
{
  return getElementBySym(sym).Z;
}

string Elements::charge2sym(int charge)
{
  return getElementByZ(charge).sym;
}

} // namespace minichem
