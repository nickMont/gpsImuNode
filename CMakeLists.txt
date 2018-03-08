cmake_minimum_required(VERSION 2.8.3)
project(gps_imu)

find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  std_msgs
  tf2_ros
  cmake_modules
  gbx_ros_bridge_msgs
)

find_package(Eigen REQUIRED)

catkin_package(
    INCLUDE_DIRS include
    CATKIN_DEPENDS gbx_ros_bridge_msgs
    DEPENDS Eigen
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
)

add_executable(gpsImu_node src/filter.cpp src/filterTW.cpp src/gpsImuNode.cpp src/gnssImuFilter.cpp src/compFilterMatrices.cpp src/gpsCallbacks.cpp)
add_dependencies(gpsImu_node gbx_ros_bridge_msgs_gencpp)
add_dependencies(gpsImu_node gbx_ros_bridge_msgs)
target_link_libraries(gpsImu_node ${catkin_LIBRARIES})

