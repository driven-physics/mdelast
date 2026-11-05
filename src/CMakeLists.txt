cmake_minimum_required(VERSION 3.0)
project(md LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(pybind11 REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED)

file(GLOB sources *.cpp *.hpp)

set(
    COMPILATION_FLAGS
    "-Wall -Wextra -Werror -pedantic-errors -Wshadow -fno-common -fstack-protector-all -rdynamic -Wpedantic -Wdouble-promotion -fstrict-aliasing"
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILATION_FLAGS}")

# Eigen emits a maybe-uninitialized warning under optimization.
set(
    CMAKE_CXX_FLAGS_RELEASE
    "${COMPILATION_FLAGS} -DNDEBUG -O2 -Wno-maybe-uninitialized -ffloat-store"
)

include_directories(${Boost_INCLUDE_DIRS})

message("pybind11 directory: " ${pybind11_DIR})

# The source filename, PYBIND11_MODULE name, and import name must match.
set(MODULE_NAME "md_driver_module")
pybind11_add_module(${MODULE_NAME} ${sources})
target_link_libraries(${MODULE_NAME} PRIVATE Eigen3::Eigen pybind11::module)
set_target_properties(
    ${MODULE_NAME}
    PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}"
)