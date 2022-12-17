# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # Alvaro Tolosa Delgado 2022
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.22)
project(BolasLocas)
# For profiling
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

set( CMAKE_BUILD_TYPE debug ) #release MinSizeRel
include_directories(include)
add_executable(main main.cpp)

add_library( tvector SHARED src/TVector.cpp)
add_library( tparticle SHARED src/TParticle.cpp)
add_library( tplot SHARED src/TPlot.cpp)
add_library( tintegrator SHARED src/TIntegrator.cpp)
target_link_libraries( tintegrator tparticle tvector)

target_link_libraries(main tvector tparticle tplot tintegrator)

set_property(TARGET main PROPERTY CXX_STANDARD 11)

find_program(GNUPLOT NAMES   gnuplot pgnuplot PATH  /bin  /usr/bin   /usr/local/bin  /sbin )

if(NOT GNUPLOT)
    message(WARNING "gnuplot not found!")
else()
    add_compile_definitions(GNUPLOT)
endif()

