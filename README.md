# C++ library for counting barcodes

![Unit tests](https://github.com/crisprVerse/kaori/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/crisprVerse/kaori/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/crisprVerse/kaori/branch/master/graph/badge.svg?token=WZkuJqiGtc)](https://codecov.io/gh/crisprVerse/kaori)

## Overview

**kaori** is a header-only C++ libary for counting the frequency of barcodes in FASTQ files.
We support a variety of barcode designs including single, combinatorial and dual barcodes in either single- or paired-end data.
Users can specify a maximum number of mismatches for identification of the target sequence (i.e., across both the constant and variable regions).
Gzipped FASTQ files can be processed, provided Zlib is available.

## Quick start

**kaori** is a header-only library, so it can be easily used by just `#include`ing the relevant source files:

```cpp
#include "kaori/kaori.hpp"
#include "byteme/byteme.hpp"

// Defining a template and the possible choices for barcodes.
std::string tmpl = "ACGT----TGCA";
std::vector<std::string> choices { "AAAA", "CCCC", "GGGG", "TTTT" };
kaori::BarcodePool pool(choices);

// Configure the handler:
// - Maximum length of 32 bp for the template.
// - Search both the forward and reverse strands.
// - Maximum mismatch of 1.
kaori::SingleBarcodeSingleEnd<32>::Options opt;
opt.max_mismatch = 1;
opt.strand = kaori::SearchStrand::BOTH;
kaori::SingleBarcodeSingleEnd<32> handler(tmpl.c_str(), tmpl.size(), pool, opt);

// Process a FASTQ file with 4 threads.
byteme::SomeFileReader input(path_to_fastq);
kaori::process_single_end(&input, handler, 4);

// Get the counts after processing.
handler.get_counts();
```

Check out the [reference documentation](https://crisprverse.github.io/kaori) for more details.

## Method description

### Algorithm

We define the "barcoding element" as the full sequence to be matched by **kaori**.
The barcoding element is parametrized as a template that consists of constant and variable regions.
Each variable region is associated with a pool of possible barcode sequences;
the template can be realized into a specific sequence of a barcoding element by replacing its variable regions with valid barcode sequences.

Given a template sequence and the associated barcode pools, **kaori** will scan each read for the template using bitwise comparisons.
If a suitable match to the template is found, the sequence of the read at each variable region is extracted and searched against the pool of known barcodes.
Imperfect barcode matches are identified through a trie-based search; the total number of mismatches is summed across the constant and variable regions.
**kaori** also caches information about imperfect matches to avoid a redundant look-up when the same sequence is encountered in later reads.

Our approach is fast and relatively easy to implement compared to full-blown sequence aligners.
Any number of mismatches are supported and the framework can be easily adapted to new barcoding configurations.
However, the downside is that indels are not supported in the search process.
We consider this limitation to be acceptable as indels are quite rare in (Illumina) sequencing data.

### Defining handlers

Each barcoding configuration is processed by a different **kaori** handler.
The code shown above for `SingleBarcodeSingleEnd` can be re-used for different handlers:

- [Single barcodes in single-end data](https://crisprverse.github.io/kaori/classkaori_1_1SingleBarcodeSingleEnd.html)
- [Single barcodes in paired-end data](https://crisprverse.github.io/kaori/classkaori_1_1SingleBarcodePairedEnd.html)
- [Combinatorial barcodes in single-end data](https://crisprverse.github.io/kaori/classkaori_1_1CombinatorialBarcodesSingleEnd.html)
- [Combinatorial barcodes in paired-end data](https://crisprverse.github.io/kaori/classkaori_1_1CombinatorialBarcodesPairedEnd.html)
- [Dual barcodes in paired-end data](https://crisprverse.github.io/kaori/classkaori_1_1DualBarcodesPairedEnd.html), 
  with [diagnostics](https://crisprverse.github.io/kaori/classkaori_1_1DualBarcodesPairedEndWithDiagnostics.html)
- [Dual barcodes in single-end data](https://crisprverse.github.io/kaori/classkaori_1_1DualBarcodesSingleEnd.html), 
  with [diagnostics](https://crisprverse.github.io/kaori/classkaori_1_1DualBarcodesSingleEndWithDiagnostics.html)
- [Random barcodes in single-end data](https://crisprverse.github.io/kaori/classkaori_1_1RandomBarcodeSingleEnd.html)

The library exports a number of utilities to easily construct a new handler - 
see the [`process_data.hpp`](https://crisprverse.github.io/kaori/process__data_8hpp.html) documentation for the handler expectations.
This can be used to quickly extend **kaori** to handle other configurations.
If you have a configuration that is not supported here, create an issue and we'll see what we can do. 
(Or even better, a pull request.)

### Setting maximum template lengths

The bitwise comparison for the constant template requires a compile-time specification of the maximum template length.
In our applications, we use templating to dispatch across a set of possible template lengths.
This improves efficiency for shorter templates while still retaining support for larger templates.
For example:

```cpp
template<size_t max_size>
std::vector<int> process(byteme::Reader* input, const std::string& tmpl, const kaori::BarcodePool& pool) {
    kaori::SingleBarcodeSingleEnd<max_size> handler(tmpl.c_str(), tmpl.size(), 0, pool, 1);
    kaori::process_single_end(input, handler, 4);
    return handler.get_counts();
}

std::vector<int> output;
if (tmpl.size() <= 32) {
    output = process<32>(&input, tmpl)
} else if (template.size() <= 64) {
    output = process<64>(&input, tmpl)
} else if (template.size() <= 128) {
    output = process<128>(&input, tmpl)
} else if (template.size() <= 256) {
    output = process<256>(&input, tmpl)
} else {
    throw std::runtime_error("template sequence is too large");
}
```

## Building projects 

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  kaori
  GIT_REPOSITORY https://github.com/crisprVerse/kaori
  GIT_TAG master # or any version of interest 
)

FetchContent_MakeAvailable(kaori)
```

Then you can link to **kaori** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe kaori)

# For libaries
target_link_libraries(mylib INTERFACE kaori)
```

### CMake with `find_package()`

You can install the library by cloning a suitable version of this repository and running the following commands:

```sh
mkdir build && cd build
cmake .. -DKAORI_TESTS=OFF
cmake --build . --target install
```

Then you can use `find_package()` as usual:

```cmake
find_package(crisprverse_kaori CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE crisprverse::kaori)
```

### Manual

If you're not using CMake, the simple approach is to just copy the files - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
You'll also need to link to the [**byteme**](https://github.com/LTLA/byteme) header-only library as well as Zlib.
