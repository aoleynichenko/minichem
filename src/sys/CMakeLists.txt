add_library(sys
    fastio.c
    time.c
)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET sys PROPERTY C_STANDARD 99)
