add_library(linalg
    linalg.c
)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET linalg PROPERTY C_STANDARD 99)

