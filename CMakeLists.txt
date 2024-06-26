cmake_minimum_required(VERSION 3.15...3.28)
project(${SKBUILD_PROJECT_NAME} 
VERSION ${SKBUILD_PROJECT_VERSION} 
LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG QUIET)

if(NOT pybind11_FOUND)
    add_subdirectory(lib/pybind11)
endif()

find_package(Matplot++ QUIET)

pybind11_add_module(projekt main.cpp)

if(NOT Matplot++_FOUND)
    add_subdirectory(lib/matplotplusplus)
    target_link_libraries(projekt PUBLIC matplot)
endif()

if(Matplot++_FOUND)
    target_link_libraries(projekt PUBLIC Matplot++::matplot)
endif()

install(TARGETS projekt LIBRARY DESTINATION .)