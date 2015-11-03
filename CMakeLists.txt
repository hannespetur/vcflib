cmake_minimum_required (VERSION 2.8)
project (vcflib)

# C++11
set (CMAKE_CXX_FLAGS "-W -Wall -Wfatal-errors -pedantic -Wno-variadic-macros -lrt -O3 -std=c++11 -DSEQAN_HAS_ZLIB=1 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/utilities")

# configure a header file to pass some of the CMake settings
# to the source code
# configure_file (
#   "${PROJECT_SOURCE_DIR}/src/constants.hpp.in"
#   "${PROJECT_BINARY_DIR}/constants.hpp"
#   )

include_directories("${PROJECT_BINARY_DIR}")

# Find external libraries
message (STATUS "Checking for zlib")
find_package(ZLIB REQUIRED)
include_directories( ${ZLIB_INCLUDE_DIRS} )

# Fetch SeqAn
find_package(SeqAn REQUIRED)

message("SeqAn include path is: ${SEQAN_INCLUDE_PATH}")
include_directories(${SEQAN_INCLUDE_DIRS})
include_directories(${SEQAN_INCLUDE_PATH})

set(EXECUTABLE_OUTPUT_PATH ../bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Find htslib
find_package(htslib REQUIRED)
include_directories(${htslib_INCLUDE_DIRS})
include_directories($ENV{HTSLIB_ROOT})

add_subdirectory(fastahack)
add_subdirectory(tabixpp)
add_subdirectory(smithwaterman)

add_subdirectory(src)

# Testing
# set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/test)
# include_directories(${EXT_PROJECTS_DIR})
# add_subdirectory(${EXT_PROJECTS_DIR}/catch)
# include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
# enable_testing(true)

# add_subdirectory(tests)
