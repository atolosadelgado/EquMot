# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # Alvaro Tolosa Delgado 2022
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.22)
project(BolasLocas)

include_directories(include)
add_executable(main main.cpp)

add_library( tvector src/TVector.cpp)
add_library( tparticle src/TParticle.cpp)
add_library( tplot src/TPlot.cpp)
add_library( tintegrator src/TIntegrator.cpp)
target_link_libraries( tintegrator tparticle tvector)

target_link_libraries(main tvector tparticle tplot tintegrator)

set_property(TARGET main PROPERTY CXX_STANDARD 11)
