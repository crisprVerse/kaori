cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG f9ed015693c45424aa2835d79fc8e58128a24c93
)

FetchContent_MakeAvailable(byteme)
