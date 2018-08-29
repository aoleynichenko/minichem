# CMake makefile for the minichem program
# 2018 A. Oleynichenko

cmake_minimum_required(VERSION 3.1)

project(minichem)
set(C_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

############################## OpenMP ##################################
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

################################ MPI ###################################

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})   # MPI header files

########################### Standard = C99 #############################

#if(CMAKE_COMPILER_IS_GNUCC)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
#endif()

############################# Compile it ! #############################

include_directories(include)

add_subdirectory(src/input)
add_subdirectory(src/sys)
add_subdirectory(src/util)
add_subdirectory(src/scf)
add_subdirectory(src/visual)
add_subdirectory(src/aoints)

add_executable(minichem.x src/main.c)

target_link_libraries(minichem.x scf aoints input util visual)
target_link_libraries(minichem.x m ${MPI_C_LIBRARIES})
target_link_libraries(minichem.x m lapack lapacke cblas)

set_property(TARGET minichem.x PROPERTY C_STANDARD 99)

