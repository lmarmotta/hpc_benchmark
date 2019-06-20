#!/usr/bin/bash

export CC=/usr/bin/mpicc
export FC=/usr/bin/mpif90
export CXX=/usr/bin/mpicxx

# HDF5 compiled with fortran and parallel support.
# FC=/usr/bin/mpif90 CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx ./configure --prefix=/home/leonardo/opt/hdf5/ --enable-parallel --enable-fortran --enable-cxx --enable-unsupported
export HDF5_DIR=/home/leonardo/opt/hdf5

# The CGNS library compiled in my home opt folder is not working with the linker.

CXX=$CXX CC=$CC FC=$FC ./configure --prefix=/home/leonardo/opt/petsc --with-scalar-type=complex --with-cc=$CC --with-fc=$FC --with-hdf5=1 --with-hdf5-dir=$HDF5_DIR --with-fortran-bindings=1 --with-cgns=1 --with-cxx=$CXX --with-blaslapack=1 --with-blas-lib=/usr/lib/libblas.so --with-lapack-lib=liblapack.so --with-valgrind=1

# Petsc will dump some instructions, so follow them ! 
# Note: Now able to install the code in a custom prefix
