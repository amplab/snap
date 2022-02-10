# SNAP

Scalable Nucleotide Alignment Program - <https://www.microsoft.com/en-us/research/project/snap/>

## Overview

SNAP is a fast and accurate aligner for short DNA reads. It is optimized for
modern read lengths of 100 bases or higher, and takes advantage of these reads
to align data quickly through a hash-based indexing scheme.

It also includes support for sorting, marking duplicates and indexing its results, eliminating the 
need for several pipeline stages used by other aligners.

## Binaries

The SNAP executable
- [v2.0.1 For Linux](https://1drv.ms/u/s!AhuEg_0yZD86hcpYCkpLlDktZnVaow?e=YX60aF)
- [v2.0.1 for Windows 10](https://1drv.ms/u/s!AhuEg_0yZD86hcpZQUgOEMrmA5qaLA?e=6FYxAv)
- [v2.0.0 for OSX](https://1drv.ms/u/s!AhuEg_0yZD86hcphrIjwoeTjdSvgoA?e=coSU85)

The SNAPCommand tool
- [SNAPCommand for Linux](https://1drv.ms/u/s!AhuEg_0yZD86hcpdvv0ZBdB1BqF57g?e=IHVbq2>)
- [SNAPCommand for Windows 10](https://1drv.ms/u/s!AhuEg_0yZD86hcpaSLKPRGJ6dcvVgA?e=vXH8y6)
- [SNAPCommand for OSX](https://1drv.ms/u/s!AhuEg_0yZD86hcpgy-ONBaw0DjFpTQ?e=cMc6eE)


## Documentation

SNAP has a one page [Quick Start Guide](https://1drv.ms/b/s!AhuEg_0yZD86hcpcvhSwRyDwk1Ru0Q?e=uAMJXV) and a more extensive [Manual](https://1drv.ms/b/s!AhuEg_0yZD86hcpblUt-muHKYsG8fA?e=R8ogug).

## Building

SNAP runs on Windows, Linux and OSX.

For Windows, we provide a Visual C++ project, `snap.sln`. Requirements:
- Visual Studio 2022

When you build it, you will have to set it to build for x64, not "Any CPU" or 32 bit.

For Linux, simply type `make`. Requirements:
- g++ version 4.8.5 or later
- zlib 1.2.11 or later from http://zlib.net/


