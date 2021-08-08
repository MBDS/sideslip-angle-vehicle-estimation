cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 14)  # Require C++14
project(fg-estimate-slip-angle)

find_package(GTSAM REQUIRED)
find_package(mrpt-math REQUIRED)
find_package(mrpt-tclap REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} mrpt::math mrpt::tclap gtsam)
