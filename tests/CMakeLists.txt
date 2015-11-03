cmake_minimum_required(VERSION 2.8)

include_directories("Including htslib2: ${htslib_INCLUDE_DIRS}")

add_executable(test_vcflib tests.cpp)

target_link_libraries(test_vcflib vcflib_core fastahack_core tabix_core smithwaterman_core ${HTSlib_LIBRARY} ${ZLIB_LIBRARIES})

target_link_libraries(test_vcflib ${LIBS})
