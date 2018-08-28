add_library(scf
    1e.c
    2e.c
    analysis.c
    boys.c
    diis.c
    eht.c
    ints_print.c
    molden.c
    rhf.c
    scf.c
    scf_input.c
    uhf.c
)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET scf PROPERTY C_STANDARD 99)
