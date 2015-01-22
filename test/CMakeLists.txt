cmake_minimum_required(VERSION 2.8.4)
project(BreakoutDetectionMain)

include_directories("${PROJECT_SOURCE_DIR}/src")
add_subdirectory("${PROJECT_SOURCE_DIR}/src")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} std=c++11 -O3")
set(CMAKE_EXE_LINKER_FLAGS "-lstdc++")
set(SOURCE_FILES test/main.cpp)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  "${PROJECT_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
add_executable(BreakoutDetectionTest ${SOURCE_FILES})
target_link_libraries(BreakoutDetectionTest breakoutdetection)