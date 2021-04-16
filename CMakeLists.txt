cmake_minimum_required(VERSION 3.11)
project(Storm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_library(PYTHIA8_LIB NAMES pythia8 HINTS "/usr/lib")
find_package(Boost REQUIRED)
find_package(fmt REQUIRED)
find_package(Threads REQUIRED)
find_package(GTest CONFIG REQUIRED)
find_package(pybind11)
find_package(GSL REQUIRED)
# find_package(Eigen3 REQUIRED)
# find_package(blaze REQUIRED)
# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fPIC")

#============================================================================
#---- Construct the library -------------------------------------------------
#============================================================================

include_directories(
        ${Boost_INCLUDE_DIRS}
        /usr/local/include
        ${CMAKE_SOURCE_DIR}/include)

set(LIB storm)

add_library(${LIB} STATIC
	# Simple model
        "${CMAKE_SOURCE_DIR}/src/simple/widths.cpp"
        "${CMAKE_SOURCE_DIR}/src/simple/msqrd.cpp"
        "${CMAKE_SOURCE_DIR}/src/simple/partially_integrated_msqrd.cpp"
        "${CMAKE_SOURCE_DIR}/src/simple/spectra.cpp"
        "${CMAKE_SOURCE_DIR}/src/pythia_spectra.cpp"
        "${CMAKE_SOURCE_DIR}/src/rambo.cpp")

target_link_libraries(${LIB} PRIVATE
        fmt::fmt
        GSL::gsl
        GSL::gslcblas
        ${Boost_LIBRARIES})

target_include_directories(${LIB} PUBLIC
        "${CMAKE_SOURCE_DIR}/include"
        ${Python3_INCLUDE_DIRS})


#============================================================================
#---- Build testing files ---------------------------------------------------
#============================================================================

set(TEST_FILES "test_model" "test_spectrum" "test_pythia")

foreach (tfile ${TEST_FILES})
    add_executable(${tfile} "tests/${tfile}.cpp")
    target_link_libraries(${tfile} PUBLIC
            ${LIB}
            ${Boost_LIBRARIES}
            Threads::Threads
            ${PYTHIA8_LIB}
            GTest::gtest_main)
endforeach ()

#============================================================================
#---- Build Python Modules ==================================================
#============================================================================

include_directories(
        ${pybind11_INCLUDE_DIRS}
        "${CMAKE_SOURCE_DIR}/include"
        ${Python3_INCLUDE_DIRS})

pybind11_add_module(_simple_widths "${CMAKE_SOURCE_DIR}/storm/models/simple/widths.cpp")
target_link_libraries(_simple_widths PRIVATE
        ${LIB}
        ${Boost_LIBRARIES}
        Threads::Threads
        ${PYTHIA8_LIB}
        ${pybind11_LIBRARIES})

pybind11_add_module(_simple_spectra "${CMAKE_SOURCE_DIR}/storm/models/simple/spectra.cpp")
target_link_libraries(_simple_spectra PRIVATE
        ${LIB}
        ${Boost_LIBRARIES}
        Threads::Threads
        ${PYTHIA8_LIB}
        ${pybind11_LIBRARIES})