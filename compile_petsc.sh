#!/usr/bin/bash

# Export basic compiling variables.

export CC=/usr/bin/mpicc
export FC=/usr/bin/mpif90
export CXX=/usr/bin/mpicxx

# HDF5 compiled with fortran and parallel support.

export HDF5_DIR=/home/leonardo/opt/hdf5

CXX=$CXX CC=$CC FC=$FC ./configure --prefix=/home/leonardo/opt/petsc --with-scalar-type=complex --with-cc=$CC --with-fc=$FC --with-hdf5=1 --with-hdf5-dir=$HDF5_DIR --with-fortran-bindings=1 --with-cgns=1 --with-cxx=$CXX --with-blaslapack=1 --with-blas-lib=/usr/lib/libblas.so --with-lapack-lib=liblapack.so --with-valgrind=1 

