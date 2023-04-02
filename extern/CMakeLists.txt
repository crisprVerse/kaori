cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

include(FetchContent)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG c638112d307eb873ec0c719046cf202d061dc4d0
)

FetchContent_MakeAvailable(byteme)
