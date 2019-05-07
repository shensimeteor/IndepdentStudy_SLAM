cmake_minimum_required(VERSION 2.8)

project(robot_pose_mle)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(Ceres REQUIRED)
find_package(Gflags REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

add_definitions(-DCERES_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})
# helloworld
add_executable(robot_pose_mle.x robot_pose_mle.cc)
target_link_libraries(robot_pose_mle.x ${CERES_LIBRARIES})
