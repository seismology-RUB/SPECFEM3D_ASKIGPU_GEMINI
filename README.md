
# SPECFEM3D_ASKIGPU_GEMINI

This is a version of [SPECFEM3D Cartesian](https://github.com/SPECFEM/specfem3d) that is able to produce
Fourier transformed wavefield output for later use in seismic full waveform inversion 
with the [ASKI software](https://github.com/seismology-RUB/ASKI_SEM3D_GEMINI).

It also accepts wavefields computed with [GEMINI](https://github.com/seismology-RUB/GEMINI_UNIFIED)
for injection at the boundaries of the SPECFEM mesh allowing hybrid forward modeling of teleseismic wavefields
through regional domains. Computation of Fourier coefficients at selected frequencies is done on the fly during
the SPECFEM time loop and can optionally be carried out on the GPU.

An detailed documentation on how to use this version of `SPECFEM` together with `ASKI` and `GEMINI` provided with the 
[ASKI software](https://github.com/seismology-RUB/ASKI_SEM3D_GEMINI). 
Python scripts helping with organising SPECFEM runs can also be found there.

This branch was initiated from a clone of SPECFEM at **commit `e4e07edc` on March 25, 2020**.

### Contributors: 
Wolfgang Friederich, Thomas Möller, Manuel Ditz, Kasper D. Fischer

## Modifications, adaptations and extensions

- To use SPECFEM 3D Cartesian for regional spherical domains, an optional mapping of the original internal Cartesian mesh to a spherical chunk was introduced in the meshing part (`create_interfaces_mesh.f90`). The Cartesian mesh is chosen as zero-centered and symmetric with respect to the $x$ and $y$-coordinates. A point $(x, y, z)$ of this Cartesian mesh is mapped to a point in a spherical chunk centered on the equator and zero longitude. Spherical coordinates of this point are $(r=R_e+z, \beta=y/R_e, \phi=x/R_e)$ with $\beta$ latitude, $\phi$ longitude, $r$ radius and $R_e$ earth radius. This choice implies that grid points with negative $z$-value are mapped to points in the earth's interior. The new Cartesian coordinate system to be used by SPECFEM is defined with origin at zero latitude and longitude and $r=R_e$. The $x$-axis points east, the $y$-axis north and the $z$-axis vertically up. The Cartesian coordinates of a point in the spherical chunk are calculated according to $(x_c=r\cos\beta\sin\phi, y_c=r\sin\beta, z_c=r\cos\beta\cos\phi-R_e)$.

- [Monteiller et al. (2013)](https://doi.org/10.1093/gji/ggs006) implemented injection of incident wavefields at all non-surface  boundary nodes of the SPECFEM mesh computed  with [AXISEM](https://github.com/geodynamics/axisem). Their code was extended here to allow injection of wavefields computed with [GEMINI](https://github.com/seismology-RUB/GEMINI_UNIFIED). To this end, coordinates and normal vectors at SPECFEM's GLL-points on the lateral and bottom surfaces of the mesh are extracted and converted to spherical coordinates. Then, for each considered seismic event, injection seismograms (particle velocity and stress vectors on the boundaries) are calculated at each point using GEMINI and stored on disk for later reading by SPECFEM. During calculation, injection seismograms are filtered to the desired frequency range, time-shifted, and optionally convolved with an event-specific moment_rate-function.

- To avoid considerable overhead time when running `SPECFEM` for computation of displacement fields from various events or Green functions from various receiver locations owing to repeated reading of the mesh and model databases as well as initializing and finalizing the message passing interface (MPI), a loop over sources was built into `SPECFEM`. Source information such as the path to the injected seismograms or receiver coordinates and single force directions for Green function computations are provided via specific input files.

- Storing of displacement fields for all events and Green functions for all receivers on disk at the high sampling rate used by SPECFEM 
requires a prohibitively vast amount of disk space. Therefore, these data are compressed by carrying out a Fourier transform at selected 
  frequencies only which is done on the fly while `SPECFEM` assembles time samples of the wavefields. This extension of `SPECFEM` was 
  already coded by [Schumacher and Friederich, (2016)](https://doi.org/10.1093/gji/ggv505) for parallel execution on CPUs and was now ported for execution on GPUs. The final output of the extended and modified version of SPECFEM are the Fourier transformed wavefields at the GLL-points of the `SPECFEM` mesh (displacements and strains) and also synthetic seismograms at all receiver positions.

------------------------------------------------------------------------------
Here follows the content of the original README file of the SPECFEM repository:

# Specfem3D

SPECFEM3D_Cartesian simulates acoustic (fluid), elastic (solid), coupled acoustic/elastic, poroelastic or seismic wave propagation in any type of conforming mesh of hexahedra (structured or not.)

It can, for instance, model seismic waves propagating in sedimentary basins or any other regional geological model following earthquakes. It can also be used for non-destructive testing or for ocean acoustics


Main "historical" authors: Dimitri Komatitsch and Jeroen Tromp
  (there are currently many more!)

## Installation

Instructions on how to install and use SPECFEM3D are
available in the

- PDF manual located in directory: [doc/USER_MANUAL](doc/USER_MANUAL)

- HTML manual (latest version): [specfem3d.readthedocs.io](http://specfem3d.readthedocs.io/)


## Development

[![Actions Status](https://github.com/geodynamics/specfem3d/workflows/CI/badge.svg)](https://github.com/geodynamics/specfem3d/actions)
[![Build Status](https://travis-ci.com/geodynamics/specfem3d.svg?branch=devel)](https://travis-ci.com/geodynamics/specfem3d)
[![codecov](https://codecov.io/gh/geodynamics/specfem3d/branch/devel/graph/badge.svg)](https://codecov.io/gh/geodynamics/specfem3d)
[![Documentation Status](https://readthedocs.org/projects/specfem3d/badge/?version=latest)](https://specfem3d.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

* Actions tests: [github actions specfem3d](https://github.com/geodynamics/specfem3d/actions)

* Travis tests: [travis-ci specfem3d](https://travis-ci.com/geodynamics/specfem3d/builds)


Development is hosted on GitHub in the
[geodynamics/specfem3d repository](https://github.com/geodynamics/specfem3d).

To contribute, please follow the guidelines located on specfem3d github wiki:
[specfem3d wiki](https://github.com/geodynamics/specfem3d/wiki)


## Computational Infrastructure for Geodynamics (CIG)

Seismology software repository: [SPECFEM3D](https://geodynamics.org/cig/software/specfem3d/)
