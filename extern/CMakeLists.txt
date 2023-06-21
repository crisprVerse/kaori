cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG a270a9cc0431d3f22e378becdd0ed1a32b7a4812
)

FetchContent_MakeAvailable(byteme)
