[![Build Status](https://travis-ci.org/Quantum-Factory/RayTransferMatrices.jl.svg?branch=master)](https://travis-ci.org/Quantum-Factory/RayTransferMatrices.jl)
[![codecov](https://codecov.io/gh/Quantum-Factory/RayTransferMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Quantum-Factory/RayTransferMatrices.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://quantum-factory.de/open-source/RayTransferMatrices.jl)

# RayTransferMatrices.jl

A Julia package for performing calculations with the [ray transfer
matrix (or ABCD)
formalism](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis),
for both 1D ray tracing and [Gaussian
beam](https://en.wikipedia.org/wiki/Gaussian_beam) propagation in the
[paraxial
approximation](https://en.wikipedia.org/wiki/Paraxial_approximation).

![](plots-02.svg)

The following introduction to the package assumes familiarity with the
ABCD formalism and its utility in optical analysis and design.  In
addition to the above links, the following are classic and useful
introductory references:

1. H. Kogelnik and T. Li, "Laser Beams and Resonators", *Applied Optics* **5**, 1550-1567 (1966)
2. A. E. Siegman, *Lasers* (University Science Books, Sausalito, 1986)

# ABCDBeamTrace.jl

This package, RayTransferMatrices.jl, is an (improved) fork of
[ABCDBeamTrace.jl](https://github.com/ngedwin98/ABCDBeamTrace.jl). It
inherits its MIT "Expat" License.

