# AssemblyCpp v5

https://arxiv.org/abs/2410.09100#

This repository contains the C++ implementation of v5 of the assembly algorithm, in a "script" style version

School of Chemistry, The University of Glasgow, University Avenue, Glasgow G12 8QQ, United Kingdom

Authors Ian Seet, Leroy Cronin


This script branch must be compiled directly from the v5/main.cpp file with the Boost Graph Library (version 1.8 or greater should be sufficient)

Compilation on both Windows and Linux should be possible with GCC or MSVC e.g.

`g++ v5/main.cpp -O AssemblyCpp -std=c++23 -o3 -I/<path_to_boost_graph_library>`

The main branch can be directly installed with cmake and does not require the Boost Graph Library

The unitTester.py script can be invoked by first navigating to the unitTests folder and calling:

`python ./unitTester.py ./<path_to_assemblyCpp_executable> batteryTest2 batteryTest2Base`
