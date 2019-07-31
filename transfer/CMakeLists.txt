cmake_minimum_required(VERSION 3.10.2)
project(MPC)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
find_package(CASADI REQUIRED)
include_directories(${CASADI_INCLUDE_DIR})
add_executable(${PROJECT_NAME} mpc_predictor_casadi.cc)
target_link_libraries(${PROJECT_NAME} ${CASADI_LIBRARIES})

