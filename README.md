![fsparse](media/logo.png)
[![DOI](https://zenodo.org/badge/686752490.svg)](https://zenodo.org/badge/latestdoi/686752490)

FSPARSE: Fortran Sparse Gallery library
=======================================

FSparse is an initiative to build a OOP API of sparse matrices with some basic kernels and utility functions using Modern Fortran.

This project is very much a work in progress, contributions are wellcome.

Available Matrices
==================
| Matrix type | no data | real | complex |
|-------------|---------|------|---------|
| COO | `COO_t` | `COO_?p` | `COO_c?p` |
| CSR | `CSR_t` | `CSR_?p` | `CSR_c?p` |
| CSC | `CSC_t` | `CSC_?p` | `CSC_c?p` |
| ELL | `ELL_t` | `ELL_?p` | `ELL_c?p` |
| SELL-C | `SELLC_t` | `SELLC_?p` | `SELLC_c?p` |

COO: COordinate Sparse format

CSR: Compressed Sparse Row format

CSC: Compressed Sparse Column format

ELL: ELLPACK

SELL-C: sliced ELLPACK format

(Where `?` stands for the precision s,d,q)

Available Kernels
==================
### Conversion
Conversion subroutines follow the naming pattern _sourcetype2targettype_m ex:
```fortran
call dense2coo( dense , coo )
```
| Matrix | dense | COO   | CSR   |
|--------|-------|-------|-------|
| dense  |       | ✅    |       |
| COO    | ✅    |       | ✅   | 
| CSR    |       | ✅    |       |
| CSC    |       |       |       |
| ELL    |       |       |       |
| SELL-C |       |       | ✅    |

### Sparse Matrix-Vector product
(available) Matrix vector products are interfaced by the procedure
```fortran
call matvec( Mat , vec_x, vec_y ) ! vec_y = vec_y + Mat * vec_x
```

| Matrix | full | symmetric |
|--------|-------|------------------|
| COO    | ✅ | ✅ |
| CSR    | ✅ | ✅ |
| CSC    | ✅ | (?) |
| ELL    | ✅ | ❌ |
| SELL-C | ✅ | ❌ |

A taste of FSPARSE
==================
```fortran
use fsparse
use iso_fortran_env, only: sp=>real32
type(COO_sp) :: COO
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
type(COO_dp) :: CSR
real(dp), allocatable :: vec_x(:)
real(dp), allocatable :: vec_y(:)

call CSR%malloc(4,5,10)
CSR%data(:)   = dble([9,-3,4,7,8,-1,8,4,5,6])
CSR%col(:)    = [1,5,1,2,2,3,4,1,3,4]
CSR%rowptr(:) = [1,3,5,8,11]

allocate( vec_x(5) , source = 1._dp )
allocate( vec_y(4) , source = 0._dp )
call matvec( CSR , vec_x, vec_y )
```
Building and using
==================
[FYPP](https://fypp.readthedocs.io/en/stable/fypp.html) is used to enable a generic programming framework, strongly inspired by the [stdlib](https://github.com/fortran-lang/stdlib) approach. A simple script is distributed to preprocess all `.fypp` files into `.f90` before building.

Only simple and double preceision
```
python deployement.py
```
Adding quad precission support
```
python deployement.py --with_qp
```

The project was built using the [Fortran Package Manager](https://github.com/fortran-lang/fpm). A manifest file is included to build and test with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `fsparse` within your FPM project, add the following to your `fpm.toml` file:
```toml
[dependencies]
fsparse = { git="https://github.com/jalvesz/FSPARSE" }
```

Inspiration & References
========================
[Iterative Methods for Sparse Linear Systems](https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf)

[Efficient Sparse Matrix-Vector Multiplication on cuda](https://www.nvidia.com/docs/io/66889/nvr-2008-004.pdf)

[gsl sparse matrices](https://www.gnu.org/software/gsl/doc/html/spmatrix.html)

[Calcul Scientifique Parallèle](https://www.dunod.com/sciences-techniques/calcul-scientifique-parallele-cours-exemples-avec-openmp-et-mpi-exercices-0)

[Implementing a Sparse Matrix Vector Product for the SELL-C/SELL-C-σ formats on NVIDIA GPUs](https://library.eecs.utk.edu/storage/files/ut-eecs-14-727.pdf)

Authors and contributors  
========================

+   [José R. Alves Z.](https://www.researchgate.net/profile/Jose-Alves-25)  
    +   Mechanical Engineer, Researcher, Scientific Software Developer
+   [Samuele Giuli](https://github.com/SamueleGiuli)
+   [Ivan Privec](https://github.com/ivan-pi)

Acknowledgement
===============

Compilation of this library was possible thanks to [Transvalor S.A.](https://www.transvalor.com/en/homepage) research activities