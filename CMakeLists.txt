cmake_minimum_required(VERSION 2.8.12)

project(kallisto)

include(GNUInstallDirs)

set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)
set(CMAKE_CXX_FLAGS_PROFILE "-g")

# Set Release type for builds where CMAKE_BUILD_TYPE is unset
# This is usually a good default as this implictly enables
#
#   CXXFLAGS = -O3 -DNDEBUG
#
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(${CMAKE_VERSION} VERSION_LESS 3.1)
    # CMake versions before 3.1 do not understand CMAKE_CXX_STANDARD
    # remove this block once CMake >=3.1 has fixated in the ecosystem
    if(APPLE)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-std=c++11>)
    else(APPLE)
        add_compile_options(-std=c++11)
    endif(APPLE)
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

#add_compile_options(-Wall -Wno-unused-function)

if(LINK MATCHES static)
    message("static build")
ELSE(LINK MATCHES shared)
    message("shared build")
ENDIF(LINK MATCHES static)

# add_compile_options(-Wdeprecated-register)

add_subdirectory(src)
include_directories(${EXT_PROJECTS_DIR})

if (BUILD_TESTING)
    add_subdirectory(${EXT_PROJECTS_DIR}/catch)

    # Includes Catch in the project:
    include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
    enable_testing(true)  # Enables unit-testing.

    add_subdirectory(unit_tests)
endif(BUILD_TESTING)

# enable_testing()
# add_test(MainTest test/tests)
