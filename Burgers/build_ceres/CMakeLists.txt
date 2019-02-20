cmake_minimum_required(VERSION 2.8)

project(Ceres_Burgers)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

FIND_PATH(CONFIG++_INCLUDE_DIR libconfig.h++ /usr/include /usr/local/include)

FIND_LIBRARY(CONFIG++_LIBRARY NAMES config++ PATH /usr/lib /usr/local/lib) 
#find_package(CONFIG++ REQUIRED)
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})



#cartrack
add_executable(Ceres_Burgers ../src/Ceres_Burgers.cpp ../src/Burgers.cpp ../src/Observation.cpp)
add_executable(run_burgers_slam.x ../src/run_burgers_slam.cpp ../src/burgers_slam.cpp ../src/Burgers.cpp ../src/Observation.cpp)
target_link_libraries(Ceres_Burgers ${CERES_LIBRARIES} ${CONFIG++_LIBRARY})
target_link_libraries(run_burgers_slam.x ${CERES_LIBRARIES} ${CONFIG++_LIBRARY})

