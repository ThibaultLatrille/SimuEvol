# SimuEvol

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/ThibaultLatrille/SimuEvol.svg?branch=master)](https://travis-ci.org/ThibaultLatrille/SimuEvol)
[![codecov](https://codecov.io/gh/ThibaultLatrille/SimuEvol/branch/master/graph/badge.svg)](https://codecov.io/gh/ThibaultLatrille/SimuEvol)

Requirements: Clang (or g++)
```bash
sudo apt install g++-5 clang-3.6
```

## Get SimuEvol up and running on Linux

### How to download and build

To get SimuEvol from a machine connected to the internet, type in a terminal:
```bash
git clone https://github.com/ThibaultLatrille/SimuEvol.git
```

This should create a folder called `SimuEvol` (the SimuEvol root folder). You must go there before compiling SimuEvol:

```bash
cd SimuEvol
```

Then, to build SimuEvol simply run:

```bash
make
```

### How to run SimuEvol

Basic usage for SimuEvol is (from the SimuEvol root folder):

```bash
bin/SimuRelax --help
bin/SimuEvol --help
bin/SimuPoly --help
```


### Authors

Thibault Latrille 