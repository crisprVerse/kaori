include(FetchContent)

FetchContent_Declare(
  byteme
  GIT_REPOSITORY https://github.com/LTLA/byteme
  GIT_TAG master # ^1.0.0
)

FetchContent_MakeAvailable(byteme)
