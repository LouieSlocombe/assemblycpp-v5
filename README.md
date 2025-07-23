# AssemblyCpp v5

This repository contains the C++ implementation of v5 of the assembly algorithm


## Installation and usage
To generate documentation, ensure doxygen is installed and run `doxygen Doxyfile` from the root. Documentation can then be accessed at `docs/html/index.html`
for installation using cmake, ensure cmake is installed and run the following from the root. This will build by default in release mode except for MSCV compiler
```
cmake -S . -B build
cmake --build build
```

When building with MSVC, you should specify release mode. 

```
cmake --build build --config release
```

This will write the assembly executable to the `build/bin` folder, except in the case of MSVC compiler, in which case it will be `build/bin/Release`

To run, e.g. on `aspirin.mol` in the same folder as the exectuable, run `./assembly aspirin` (note the .mol is ommitted)

The assembly index is output to stdout and to the file aspirinOut. Full pathway details are output to the file aspirinPathway. If terminated with Ctrl-C, the program will output the best assembly index found so far, provided the enumeration process is complete.

Cmake will also attempt to install the catch2 testing suite and compile the `unit_tests` binary in `build/bin`. This currently includes integration tests run on the assembly binary (in its compilation location of `build/bin`), and a couple of trivial unit tests, with more unit tests to be added. 

The integration test checks completion and expected assembly index for ~1k molecules. Failure details are output to integration_failures.txt in the root (only generated if there are failures). A file integration_timings.csv will also be generated in the root, which includes timings that can be used for basic benchmarking, although many of the molecules run very quickly so the utility may be limited.

A further speed test can be run using a python script in `tests/speed/speedtest.py`. The location of the assembly executable should be passed to the script, and the file should be run from the `tests/speed` folder. If the assembly executable is in the default `build/bin` then run

```
python speedtest.py ../../build/bin/assembly
```

or on Windows

```
python speedtest.py ../../build/bin/assembly.exe
```