#!/usr/bin/env python

# Test suite for the minichem package

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from minitest import Test, Filter

########################################################################
#                            S     C     F
########################################################################

print '>>> RHF'

# CO
t1_scf  = Filter("Total SCF energy = ", -111.22544951, 1e-7)
t1_enuc = Filter("Nuclear repulsion energy = ", 22.17458792, 1e-7)
t1_dip  = Filter("|D| = ", 0.04893493, 1e-7)
Test("CO molecule RHF/STO-3G, symm C1", "CO.nw", filters=[t1_scf,t1_enuc,t1_dip]).run()

# H2O
t1_scf  = Filter("Total SCF energy = ", -74.96314638, 1e-7)
t1_enuc = Filter("Nuclear repulsion energy = ", 9.18922278, 1e-7)
t1_dip  = Filter("|D| = ", 0.67964434, 1e-6)
Test("H2O molecule dir-RHF/STO-3G, symm C1", "H2O.nw", filters=[t1_scf,t1_enuc,t1_dip]).run()

# C6H6
t1_scf  = Filter("Total SCF energy = ", -227.89074473, 1e-7)
t1_enuc = Filter("Nuclear repulsion energy = ", 203.36796552, 1e-7)
t1_dip  = Filter("|D| = ", 0.0, 1e-6)
Test("C6H6 molecule RHF/STO-3G, symm C1", "C6H6.nw", filters=[t1_scf,t1_enuc,t1_dip]).run()


print '>>> UHF'

# Li atom
t1_scf  = Filter("Total SCF energy = ", -7.31552589, 1e-7)
t1_enuc = Filter("Nuclear repulsion energy = ", 0.0, 1e-7)
t1_s2   = Filter("UHF <S^2> = ", 0.7500, 1e-5)
Test("Li atom UHF/STO-3G, symm C1", "Li.nw", filters=[t1_scf,t1_enuc,t1_s2]).run()

#CH3
t1_scf  = Filter("Total SCF energy = ", -39.07664260, 1e-7)
t1_enuc = Filter("Nuclear repulsion energy = ", 9.73109626, 1e-7)
t1_s2   = Filter("UHF <S^2> = ", 0.7648, 1e-5)
Test("CH3 radical UHF/STO-3G, symm C1", "CH3.nw", filters=[t1_scf,t1_enuc,t1_s2]).run()
