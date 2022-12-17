# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # Alvaro Tolosa Delgado 2022
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.22)
project(BolasLocas)
# For profiling
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

find_package(OpenMP REQUIRED COMPONENTS CXX)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


set( CMAKE_BUILD_TYPE RELEASE ) #release MinSizeRel debug
include_directories(include)
add_executable(main main.cpp)

add_library( tvector SHARED src/TVector.cpp)
add_library( tparticle SHARED src/TParticle.cpp)
add_library( tplot SHARED src/TPlot.cpp)
add_library( tmatrixtriang SHARED src/TMatrixtriang.cpp)
add_library( tintegrator SHARED src/TIntegrator.cpp)
target_link_libraries( tintegrator tparticle tvector tmatrixtriang OpenMP::OpenMP_CXX)

target_link_libraries(main tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX)


set_property(TARGET main PROPERTY CXX_STANDARD 11)

find_program(GNUPLOT NAMES   gnuplot pgnuplot PATH  /bin  /usr/bin   /usr/local/bin  /sbin )

if(NOT GNUPLOT)
    message(WARNING "gnuplot not found!")
else()
    add_compile_definitions(GNUPLOT)
endif()

