add_library(sys
    fastio.c
)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET sys PROPERTY C_STANDARD 99)
