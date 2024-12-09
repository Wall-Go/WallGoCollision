cmake_minimum_required(VERSION 3.18)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(${CMAKE_SOURCE_DIR}/cmake/GenerateProjectVersion.cmake)
GenerateProjectVersion(WALLGOCOLLISION_VERSION)

project(Collision VERSION ${WALLGOCOLLISION_VERSION} LANGUAGES CXX)

option(USE_CXX20 "Prefer C++20 if available (may allow small optimizations)" ON)

if (USE_CXX20)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-std=c++20" COMPILER_SUPPORTS_CXX20)
    
    if (COMPILER_SUPPORTS_CXX20)
        set(CMAKE_CXX_STANDARD 20)
    else()
        # Fallback to C++17
        set(CMAKE_CXX_STANDARD 17)
    endif()
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

# What to call the C++ library
set(WALLGO_LIB WallGoCollision CACHE STRING "Name of the WallGo Collision library")

set(INSTALL_DIR ${CMAKE_SOURCE_DIR}/bin)

# Set default built type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

if (CMAKE_TOOLCHAIN_FILE)
    message(STATUS "Using toolchain file: ${CMAKE_TOOLCHAIN_FILE}")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

## Compile with profiler support?
option(PROFILER "Compile with profiler flags" OFF)
if (PROFILER)
    add_compile_options(-pg)
    add_link_options(-pg)
endif()

## Use OpenMP by default
option(USE_OMP "Build with OpenMP support" ON)
option(REQUIRE_OMP "Require OpenMP and abort the build if not present" OFF)
option(BUILD_PYTHON_MODULE "Build Python module" ON)
option(BUILD_EXAMPLES "Build WallGo/Collision C++ examples" ON)

if (REQUIRE_OMP)
    set(USE_OMP ON)
endif()

if(BUILD_PYTHON_MODULE)
	add_compile_definitions(WITH_PYTHON=1)
else()
    add_compile_definitions(WITH_PYTHON=0)
endif()

# Compiler-specific settings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

    if (NOT USE_OMP)
        # Disable pragma warnings, otherwise each #pragma omp would throw a warning
        add_compile_options(-Wno-unknown-pragmas)
    endif()

    add_compile_options(-Wall -Wextra)
    #add_compile_options(-fvisibility=hidden)

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		message(STATUS "==== NOTE: Configuring Debug build")
		message(STATUS "")

        # Treat all warnings as errors?
        #add_compile_options(-Werror)
        # Allow compiler to become annoying
        add_compile_options(-Wpedantic)
    endif()
endif()

if (MSVC AND NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(/W4)
endif()

add_subdirectory(src)

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if (BUILD_PYTHON_MODULE)
	add_subdirectory(python)
endif()

## Tests are up-to-date for the python module only
#if(BUILD_UNIT_TESTS)
#    add_subdirectory(tests)
#endif()