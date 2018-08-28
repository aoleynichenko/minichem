add_library(input
    basis.c
    chem.c
    inp.c
    lexer.c
)
target_link_libraries(input scf)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET input PROPERTY C_STANDARD 99)
