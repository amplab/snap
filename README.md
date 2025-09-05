# SNAP

Scalable Nucleotide Alignment Program - <https://www.microsoft.com/en-us/research/project/snap/>

## Overview

SNAP is a fast and accurate aligner for short DNA reads. It is optimized for
modern read lengths of 100 bases or higher, and takes advantage of these reads
to align data quickly through a hash-based indexing scheme.

It also includes support for sorting, marking duplicates and indexing its results, eliminating the 
need for several pipeline stages used by other aligners.

## Binaries

Current binaries are available on this GitHub page under "releases"


## Documentation

SNAP has a one page [Quick Start Guide](https://1drv.ms/b/s!AhuEg_0yZD86hcpcvhSwRyDwk1Ru0Q?e=4BvzLn) and a more extensive [Manual](https://1drv.ms/b/s!AhuEg_0yZD86hcpblUt-muHKYsG8fA?e=mbyUP5).

## Building

SNAP runs on Windows, Linux and OSX.

For Windows, we provide a Visual C++ solution, `snap.sln`. Requirements:
- Visual Studio 2022

When you build it, you will have to set it to build for x64, not "Any CPU" or 32 bit.

For Linux, simply type `make`. Requirements:
- g++ version 4.8.5 or later
- zlib 1.2.11 or later from http://zlib.net/.  On Ubuntu, do sudo apt install libz1g-dev


