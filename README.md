# SNAP

Scalable Nucleotide Alignment Program - <http://snap.cs.berkeley.edu>

## Overview

SNAP is a fast and accurate aligner for short DNA reads. It is optimized for
modern read lengths of 100 bases or higher, and takes advantage of these reads
to align data quickly through a hash-based indexing scheme.

## Documentation

A quick start guide and user manual are available in the `docs` folder, with
additional documentation at <http://snap.cs.berkeley.edu>.

## Building

SNAP runs on Windows, Linux and Mac OS X.

For Windows, we provide a Visual C++ project, `snap.sln`. Requirements:
- Visual Studio 2012 (11.0)

For Linux and OS X, simply type `make`. Requirements:
- g++ version 4.6
- zlib 1.2.8 from http://zlib.net/


