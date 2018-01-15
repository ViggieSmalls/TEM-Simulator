cmake_minimum_required(VERSION 3.8)
project(TEM-Simulator)

include(ExternalProject)

set(FFTW3_URL http://www.fftw.org/fftw-3.3.7.tar.gz)
set(FFTW3_PREFIX fftw-3.3.7)

ExternalProject_Add(libfftw3
        URL ${FFTW3_URL}
        URL_HASH MD5=0d5915d7d39b3253c1cc05030d79ac47
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${FFTW3_PREFIX}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${FFTW3_PREFIX}
)

ExternalProject_Get_Property(libfftw3 install_dir)
include_directories(${install_dir}/include)
link_directories(${install_dir}/lib)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -pedantic -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall -pedantic -O3")

AUX_SOURCE_DIRECTORY(src SOURCE_FILES)

#include_directories(C:/mingw-w64/libs/include)
#link_directories(C:/mingw-w64/libs/lib)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES})
add_dependencies(${CMAKE_PROJECT_NAME} libfftw3)
target_link_libraries(${CMAKE_PROJECT_NAME} fftw3 m)
