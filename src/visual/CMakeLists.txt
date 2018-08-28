add_library(visual
    molden.c
)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET visual PROPERTY C_STANDARD 99)

