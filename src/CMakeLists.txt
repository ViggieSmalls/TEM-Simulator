cmake_minimum_required(VERSION 3.8)
project(TEM-Simulator)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -O3")

AUX_SOURCE_DIRECTORY(. SOURCE_FILES)

include_directories(C:/mingw-w64/libs/include)
link_directories(C:/mingw-w64/libs/lib)

add_executable(TEM-Simulator ${SOURCE_FILES})

target_link_libraries(TEM-Simulator libfftw3-3)
