cmake_minimum_required(VERSION 2.8)

# Create a vcflib
set(vcflib_SOURCE_FILES
   Variant.cpp
   split.cpp
   ssw_cpp.cpp
)

add_library(vcflib_core ${vcflib_SOURCE_FILES} )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(vcfcheck vcfcheck.cpp)
target_link_libraries(vcfcheck vcflib_core fastahack_core tabix_core smithwaterman_core ${HTSlib_LIBRARY} ${ZLIB_LIBRARIES})

# Add the executable
target_link_libraries(vcfcheck ${LIBS})

# Installs vcflib binaries to /usr/local/bin using [sudo] make install
install (TARGETS vcfcheck DESTINATION bin)
