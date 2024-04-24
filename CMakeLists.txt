cmake_minimum_required(VERSION 3.15...3.28)
project(${SKBUILD_PROJECT_NAME} 
VERSION ${SKBUILD_PROJECT_VERSION} 
LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(Matplot++ REQUIRED)

pybind11_add_module(projekt main.cpp)
target_link_libraries(projekt PUBLIC Matplot++::matplot)
#set_property(TARGET matplot PROPERTY POSITION_INDEPENDENT_CODE ON)

install(TARGETS projekt LIBRARY DESTINATION .)