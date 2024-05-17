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

add_subdirectory(lib/audiofile)

add_definitions (-DPROJECT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
file (COPY test-audio.wav DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(projekt PUBLIC AudioFile)

install(TARGETS projekt LIBRARY DESTINATION .)
