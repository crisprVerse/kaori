cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG 964743bf8f3cfcabf6fb5435c998062c5765ee3c
)

FetchContent_MakeAvailable(byteme)
