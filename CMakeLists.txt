cmake_minimum_required(VERSION 3.14)

project(kaori
    VERSION 1.0.0
    DESCRIPTION "Header-only C++ library for screen counting"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(kaori INTERFACE)

target_include_directories(kaori INTERFACE include/)

add_subdirectory(extern)

target_link_libraries(kaori INTERFACE byteme)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
