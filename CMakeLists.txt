cmake_minimum_required(VERSION 3.22)
project(BolasLocas)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(main main.cpp)
target_link_libraries(main ${ARMADILLO_LIBRARIES})
set_property(TARGET main PROPERTY CXX_STANDARD 11)
