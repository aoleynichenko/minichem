add_library(scf
    analysis.c
    diis.c
    eht.c
    rhf.c
    scf.c
    scf_input.c
    uhf.c
)

target_link_libraries(scf aoints visual linalg util)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET scf PROPERTY C_STANDARD 99)
