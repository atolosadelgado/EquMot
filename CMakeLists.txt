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
target_link_libraries( tintegrator tparticle tvector tmatrixtriang tplot OpenMP::OpenMP_CXX TBB::tbb) #

target_link_libraries(main tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb ) #

find_package(TBB REQUIRED)
target_compile_features(main PRIVATE cxx_std_20)
target_compile_features(tintegrator PRIVATE cxx_std_20)


set_property(TARGET main PROPERTY CXX_STANDARD 20)

find_program(GNUPLOT NAMES   gnuplot pgnuplot PATH  /bin  /usr/bin   /usr/local/bin  /sbin )

if(NOT GNUPLOT)
    message(WARNING "gnuplot not found!")
else()
    add_compile_definitions(GNUPLOT)
endif()

enable_testing()

add_executable(test_force_seq src/test_force_seq.cpp)
target_link_libraries(test_force_seq tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_force_seq COMMAND  $<TARGET_FILE:test_force_seq> )


add_executable(test_force_openmp src/test_force_openmp.cpp)
target_link_libraries(test_force_openmp tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_force_openmp COMMAND  $<TARGET_FILE:test_force_openmp> )



add_executable(test_force_triang_seq src/test_force_triang_seq.cpp)
target_link_libraries(test_force_triang_seq tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_force_triang_seq COMMAND  $<TARGET_FILE:test_force_triang_seq> )


add_executable(test_force_triang_openmp src/test_force_triang_openmp.cpp)
target_link_libraries(test_force_triang_openmp tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_force_triang_openmp COMMAND  $<TARGET_FILE:test_force_triang_openmp> )


add_executable(test_force_cxx src/test_force_cxx.cpp)
target_link_libraries(test_force_cxx tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_force_cxx COMMAND  $<TARGET_FILE:test_force_cxx> )

add_executable(test_force_integrator_cxx src/test_force_integrator_cxx.cpp)
target_link_libraries(test_force_integrator_cxx tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_force_integrator_cxx COMMAND  $<TARGET_FILE:test_force_integrator_cxx> )

add_executable(test_full_cxx src/test_full_cxx.cpp)
target_link_libraries(test_full_cxx tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_full_cxx COMMAND  $<TARGET_FILE:test_full_cxx> )

add_executable(test_full_tbb src/test_full_tbb.cpp)
target_link_libraries(test_full_tbb tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_full_tbb COMMAND  $<TARGET_FILE:test_full_tbb> )

add_executable(test_SolarSystem src/test_SolarSystem.cpp)
target_link_libraries(test_SolarSystem tvector tparticle tplot tintegrator tmatrixtriang OpenMP::OpenMP_CXX TBB::tbb) #   TBB::tbb
add_test( NAME test_SolarSystem COMMAND  $<TARGET_FILE:test_SolarSystem> )
