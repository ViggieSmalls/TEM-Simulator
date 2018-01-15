cmake_minimum_required(VERSION 3.5)
project(TEM-Simulator)

include(ExternalProject)

set(FFTW3_URL http://www.fftw.org/fftw-3.3.7.tar.gz)
set(FFTW3_PREFIX fftw-3.3.7)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
	set(SHARED OFF)
endif( CMAKE_BUILD_TYPE STREQUAL "Release" )

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(SHARED ON)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

ExternalProject_Add(libfftw3
        URL ${FFTW3_URL}
        URL_HASH MD5=0d5915d7d39b3253c1cc05030d79ac47
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${FFTW3_PREFIX}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${FFTW3_PREFIX} -DBUILD_SHARED_LIBS=${SHARED}
)

ExternalProject_Get_Property(libfftw3 INSTALL_DIR)
include_directories(${INSTALL_DIR}/include)
link_directories(${INSTALL_DIR}/lib)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -pedantic -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall -pedantic -O3 -static")

aux_source_directory(src SOURCE_FILES)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCE_FILES})

add_dependencies(${CMAKE_PROJECT_NAME} libfftw3)

target_link_libraries(${CMAKE_PROJECT_NAME} fftw3 m)
