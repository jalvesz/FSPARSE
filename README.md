![fsparse](media/logo.png)

FSPARSE: Fortran Sparse Gallery library
=======================================

FSparse is an initiative to build a OOP API of sparse matrices with some basic kernels and utility functions using Modern Fortran.

This project is very much a work in progress, contributions are wellcome.

Available Matrices
==================
COO: COordinate Sparse format

CSR: Compressed Sparse Row format

CSC: Compressed Sparse Column format

ELL: ELLPACK

All matrices come currently in three flavors:

type(COO_t) : no data buffer (only index targeting graph manipulation)

type(COOr32_t) : simple precision data buffer

type(COOr64_t) : double precision data buffer

A taste of FSPARSE
==================
```fortran
use fsparse
use iso_fortran_env, only: sp=>real32
type(COOr32_t) :: COO
real(sp), allocatable :: dense(:,:)

allocate( dense(4,5) )
dense = reshape([9.0,4.0,0.0,4.0,0.0, &
                 7.0,8.0,0.0,0.0,0.0, &
                 -1.0,5.0,0.0,0.0,8.0,&
                 6.0,-3.0,0.0,0.0,0.0],[4,5])

call dense2coo( dense , COO )
```

```fortran
use fsparse
use iso_fortran_env, only: dp=>real64
type(CSRr64_t) :: CSR
real(dp), allocatable :: vec_x(:)
real(dp), allocatable :: vec_y(:)

call CSR%malloc(4,5,10)
CSR%data(:)   = dble([9,-3,4,7,8,-1,8,4,5,6])
CSR%col(:)    = [1,5,1,2,2,3,4,1,3,4]
CSR%rowptr(:) = [1,3,5,8,11]

allocate( vec_x(5) , source = 1._dp )
allocate( vec_y(4) , source = 0._dp )
call matvec( CSR , vec_x , vec_y )
```

Inspiration
===========
[Efficient Sparse Matrix-Vector Multiplication on cuda](https://www.nvidia.com/docs/io/66889/nvr-2008-004.pdf)

[gsl sparse matrices](https://www.gnu.org/software/gsl/doc/html/spmatrix.html)

[Calcul Scientifique Parallèle](https://www.dunod.com/sciences-techniques/calcul-scientifique-parallele-cours-exemples-avec-openmp-et-mpi-exercices-0)

Authors and contributors  
========================

+   [José R. Alves Z.](https://www.researchgate.net/profile/Jose-Alves-25)  
    +   Mechanical Engineer, Researcher, Scientific Software Developer

Acknowledgement
===============

Compilation of this library was possible thanks to [Transvalor S.A.](https://www.transvalor.com/en/homepage) research activities
