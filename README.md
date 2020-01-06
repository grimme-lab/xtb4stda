# sTDA-xTB for ground state calculations

[![Build Status](https://travis-ci.com/grimme-lab/xtb4stda.svg?branch=master)](https://travis-ci.com/grimme-lab/xtb4stda)
[![Build Status](https://github.com/grimme-lab/xtb4stda/workflows/CI/badge.svg)](https://github.com/grimme-lab/xtb4stda/actions)

This project provides `xtb4stda`, a program to calculate the ground state with
sTDA-xTB to be used in further [`stda`](https://github.com/grimme-lab/stda)
calculations.

## Installation

Statically linked binaries can be found at the projects
[release page](https://github.com/grimme-lab/xtb4stda/releases/latest).
To build from source this project uses a `make` based build system and requires
a version of Intel Parallel Studio 17 or newer to be compiled (also requires
`ruby` for some hacky scripts).
To trigger the build run in the root directory

```bash
make
```

You will find a statically linked executable in `exe/xtb4stda`.
To make `xtb4stda` accessible export

```bash
export XTB4STDAHOME=$PWD
export PATH=$PATH:XTB4STDAHOME/exe
```

For parallel usage set the threads for OMP and the MKL linear algebra backend by

```bash
export OMP_NUM_THREADS=<ncores> MKL_NUM_THREADS=<ncores>
```

For larger systems please adjust the stack size accordingly, otherwise
stack overflows *will* occur. Use something along the lines of this:

```bash
ulimit -s unlimited
export OMP_STACKSIZE=4G
```

To make adjustments to the build system check the directory `MAKE/`.

### Alternatives

If you are not a fan of `make`, you can use [`meson`](https://mesonbuild.com/)
as alternative, but it requires a fairly new version like 0.49 or newer for a
decent Fortran support.
For the default backend [`ninja`](https://ninja-build.org/) version 1.5 or newer
has to be provided.

To perform a build run (with GCC, run `export FC=ifort` for Intel builds)

```bash
meson setup build_gcc
ninja -C build_gcc
```

This also allows to install the program locally by running

```bash
[sudo] ninja -C build_gcc install
```

Which will default to an installation in `/usr/local` and should automatically
include `xtb4stda` in your systems `PATH` variable.
The environment variable to set is than usually only

```bash
export XTB4STDAHOME=/usr/local/share/xtb4stda
```

## Usage

The `xtb4stda` binary will read the geometry from its *first* command line argument
and the input is assumed to be either a Turbomole coordinate data group or a
xyz-file. The extension of the file is not used to distinguish files.

```bash
xtb4stda coord > gs.stda-xtb.out
```

After the run you will find a `wfn.xtb` file in the directory which can be used
with the [`stda`](https://github.com/grimme-lab/stda) program.

## Citations

- S. Grimme and C. Bannwarth, *J. Chem. Phys.*, **2016**, 145, 054103.
  DOI: [10.1063/1.4959605](https://dx.doi.org/10.1063/1.4959605)

## License

`xtb4stda` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`xtb4stda` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

## Bugs

A bug is a *demonstratable problem* caused by the code in this repository.
Good bug reports are extremely valuable for us - thank you!

Before opening a bug report:

1. Check if the issue has already been reported.
2. Check if it still is an issue or has already been fixed?
   Try to reproduce it with the latest version from the `master` branch.
3. Isolate the problem and create a reduced test case.

A good bug report should not leave others needing to chase you up for more
information. So please try to be as detailed as possible in your report,
answer at least these questions:

1. Which version of `xtb4stda` are you using? The current version is always
   a subject to change, so be more specific.
2. What is your environment (your laptop, the cluster of the university)?
3. What steps will reproduce the issue?
   We have to reproduce the issue, so we need all the input files.
4. What would be the expected outcome?
5. What did you see instead?

All these details will help people to fix any potential bugs.

### Known Issues

For large systems with more than 33000 basis functions an integer overflow
in the linear algebra backend will occur. To amend this issue the integer
precision range must be increased from 32 bits to 64 bits by recompiling
the program.
Adjust the make build by changing `-lmkl_intel_lp64` to `-lmkl_intel_il64`
in [`MAKE/Makerules`](https://github.com/grimme-lab/xtb4stda/blob/master/MAKE/Makerules)
or when using the `meson` backend by adding `-Dinterface=64` in the
configuration step.
Note that this option is currently only supported with the MKL backend.
