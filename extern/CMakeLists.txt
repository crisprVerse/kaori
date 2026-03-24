include(FetchContent)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG master # ^2.1.1
)

FetchContent_MakeAvailable(byteme)
