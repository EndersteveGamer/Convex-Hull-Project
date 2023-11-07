cmake_minimum_required(VERSION 3.22)
project(Hull C)

set(CMAKE_C_STANDARD 99)

add_executable(Hull hull.c)

target_link_libraries(Hull PRIVATE m)