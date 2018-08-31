add_library(util
    decor.c
    errors.c
    memory.c
    misc.c
    timer.c
)

target_link_libraries(util sys)

set(C_STANDARD_REQUIRED ON)
set_property(TARGET util PROPERTY C_STANDARD 99)
