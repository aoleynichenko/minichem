add_library(aoints
    1e.c
    2e.c
    boys.c
    ints_print.c
)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET scf PROPERTY C_STANDARD 99)
