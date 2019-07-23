cmake_minimum_required (VERSION 3.15)

project(EKFc C)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "-O0")

add_executable(EKFc src/main.c src/main.h src/matrix.c src/matrix.h)