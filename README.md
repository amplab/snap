# SNAP

Scalable Nucleotide Alignment Program - <http://snap.cs.berkeley.edu>

## Overview

SNAP is a fast and accurate aligner for short DNA reads. It is optimized for
modern read lengths of 100 bases or higher, and takes advantage of these reads
to align data quickly through a hash-based indexing scheme.

## Binaries

[v1.0.0 For Linux](https://1drv.ms/u/s!AhuEg_0yZD86hcpYCkpLlDktZnVaow?e=YX60aF)

[v1.0.0 for Windows 10](https://1drv.ms/u/s!AhuEg_0yZD86hcpZQUgOEMrmA5qaLA?e=6FYxAv)


[SnapCommand for Linux](https://1drv.ms/u/s!AhuEg_0yZD86hcpaSLKPRGJ6dcvVgA?e=swOlYD)

[SnapCommand for Windows 10](https://1drv.ms/u/s!AhuEg_0yZD86hcpaSLKPRGJ6dcvVgA?e=vXH8y6)


## Documentation

SNAP has a one page [Quick Start Guide](https://1drv.ms/b/s!AhuEg_0yZD86hcpcvhSwRyDwk1Ru0Q?e=uAMJXV) and a more extensive [Manual](https://1drv.ms/b/s!AhuEg_0yZD86hcpblUt-muHKYsG8fA?e=R8ogug)

## Building

SNAP runs on Windows and Linux.

For Windows, we provide a Visual C++ project, `snap.sln`. Requirements:
- Visual Studio 2019

When you build it, you will have to set it to build for x64, not "Any CPU" or 32 bit.

For Linux and OS X, simply type `make`. Requirements:
- g++ version 4.6
- zlib 1.2.8 from http://zlib.net/


