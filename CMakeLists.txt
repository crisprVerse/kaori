cmake_minimum_required(VERSION 3.24)

project(kaori
    VERSION 1.1.2
    DESCRIPTION "Header-only C++ library for screen counting"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(kaori INTERFACE)
add_library(crisprverse::kaori ALIAS kaori)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

target_include_directories(kaori INTERFACE 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/crisprverse_kaori>")

option(KAORI_FETCH_EXTERN "Automatically fetch kaori's external dependencies." ON)
if(KAORI_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(ltla_byteme 1.0.0 CONFIG REQUIRED)
endif()

target_link_libraries(kaori INTERFACE ltla::byteme)

# Building the test-related machinery, if we are compiling this library directly.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(KAORI_TESTS "Build kaori's test suite." ON)
else()
    option(KAORI_TESTS "Build kaori's test suite." OFF)
endif()

if(KAORI_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()

# Installing for find_package.
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/crisprverse_kaori)

install(TARGETS kaori
    EXPORT kaoriTargets)

install(EXPORT kaoriTargets
    FILE crisprverse_kaoriTargets.cmake
    NAMESPACE crisprverse::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/crisprverse_kaori)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/crisprverse_kaoriConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/crisprverse_kaori)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/crisprverse_kaoriConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/crisprverse_kaoriConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/crisprverse_kaoriConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/crisprverse_kaori)
