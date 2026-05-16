# HEALPix Full-Sky FITS Map Analysis Toolkit

**Warning:** This repository is currently **under active development**.  
Features, APIs, file structures, and program behavior may change without notice.  
Some routines may still be experimental, incomplete, or insufficiently tested.

A collection of programs for analyzing full-sky maps in FITS format using the HEALPix pixelization scheme. This repository includes tools for component separation, signal decomposition, and general map manipulation commonly used in astrophysics and cosmology.

## Features

Current programs included in the repository:

- **Global Internal Linear Combination (ILC)**  
  Performs component separation using global ILC weights.

- **Local Internal Linear Combination (Local ILC)**  
  Performs spatially localized ILC component separation.

- **Empirical Mode Decomposition (EMD)**  
  Decomposes HEALPix maps into intrinsic mode functions (IMFs).

- **Multivariate Empirical Mode Decomposition (MEMD)**  
  Multichannel extension of EMD for simultaneous decomposition of multiple maps.

- **Map Summation Utilities**  
  Tools for combining multiple HEALPix maps.

- Additional utilities for processing and manipulating full-sky FITS maps.

## Supported Formats

- FITS full-sky maps
- HEALPix pixelization scheme

## Applications

These tools are intended for applications such as:

- Cosmic Microwave Background (CMB) component separation
- Astrophysical foreground analysis
- Multi-frequency sky map analysis
- Signal decomposition and denoising
- Statistical analysis of spherical data

## Requirements

Typical dependencies include:

- HEALPix
- CFITSIO
- Fortran compiler
- Numerical libraries (BLAS/LAPACK)
