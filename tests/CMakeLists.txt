include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Avoid installing GoogleTest when installing this project.
option(INSTALL_GTEST "Enable installation of googletest." OFF)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
    libtest 
    src/FastqReader.cpp
    src/ScanTemplate.cpp
    src/MismatchTrie.cpp
    src/BarcodeSearch.cpp
    src/SimpleSingleMatch.cpp
    src/process_data.cpp
    src/handlers/SingleBarcodeSingleEnd.cpp
    src/handlers/SingleBarcodePairedEnd.cpp
    src/handlers/CombinatorialBarcodesSingleEnd.cpp
    src/handlers/CombinatorialBarcodesPairedEnd.cpp
    src/handlers/DualBarcodesPairedEnd.cpp
    src/handlers/DualBarcodesPairedEndWithDiagnostics.cpp
    src/handlers/DualBarcodesSingleEnd.cpp
    src/handlers/DualBarcodesSingleEndWithDiagnostics.cpp
    src/handlers/RandomBarcodeSingleEnd.cpp
)

target_link_libraries(
    libtest
    gtest_main
    kaori
)

target_compile_options(libtest PRIVATE -Wall -Wextra -Wpedantic -Werror)

set(CODE_COVERAGE "Enable coverage testing" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

include(GoogleTest)
gtest_discover_tests(libtest)
