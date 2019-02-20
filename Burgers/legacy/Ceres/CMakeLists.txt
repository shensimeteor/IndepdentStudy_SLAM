cmake_minimum_required(VERSION 2.8)

project(Ceres_Burgers)

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

#cartrack
add_executable(Ceres_Burgers Ceres_Burgers.cpp Burgers.cpp Observation.cpp)
target_link_libraries(Ceres_Burgers ${CERES_LIBRARIES})


