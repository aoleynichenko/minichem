add_library(aoints
    aoints.c
    1e.c
    2e.c
    boys.c
    ints_print.c
)

target_link_libraries(aoints sys)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET scf PROPERTY C_STANDARD 99)
