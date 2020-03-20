cmake_minimum_required(VERSION 3.14)
project(crisprSDS)

set(CMAKE_CXX_STANDARD 14)
set(SOURCE_FILES main.cpp)

add_executable(crisprSDS main.cpp src/crispr.cpp includes/crispr.h)
add_executable(main ${SOURCE_FILES})

include_directories(includes)
include_directories(external/sdsl-lite/include)
include_directories(external/sdsl-lite/build/external/libdivsufsort/include)

target_link_libraries(crisprSDS ${CMAKE_SOURCE_DIR}/external/sdsl-lite/build/lib/libsdsl.a)
target_link_libraries(crisprSDS ${CMAKE_SOURCE_DIR}/external/sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort.a)
target_link_libraries(crisprSDS ${CMAKE_SOURCE_DIR}/external/sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a)


include_directories(includes)

add_subdirectory(src)