cmake_minimum_required(VERSION 2.8.4)
include_directories("${PROJECT_SOURCE_DIR}/src")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fno-strict-aliasing -fno-common -dynamic -arch i386 -arch x86_64 -DNDEBUG -fwrapv -Wall -Wstrict-prototypes")
set(SOURCE_FILES edm-per.cpp edm-multi.cpp edmTail.cpp edmx.cpp helper.cpp)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  "${PROJECT_SOURCE_DIR}/lib")
add_library(breakoutdetection STATIC ${SOURCE_FILES})