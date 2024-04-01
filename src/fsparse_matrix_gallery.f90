!---------------------------------------------------
! Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module fsparse_matrix_gallery

    !! This module provides the base sparse matrix types

    use iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private
  
    public :: COO_t, COOr32_t, COOr64_t, COOc32_t, COOc64_t
    public :: CSR_t, CSRr32_t, CSRr64_t, CSRc32_t, CSRc64_t
    public :: CSC_t, CSCr32_t, CSCr64_t, CSCc32_t, CSCc64_t
    public :: ELL_t, ELLr32_t, ELLr64_t, ELLc32_t, ELLc64_t

    
    ! -- Global parameters
    integer, parameter, public :: k_NOSYMMETRY = 0 !! Full Sparse matrix (no symmetry considerations)
    integer, parameter, public :: k_SYMTRIINF = 1  !! Symmetric Sparse matrix with triangular inferior storage
    integer, parameter, public :: k_SYMTRISUP = 2  !! Symmetric Sparse matrix with triangular supperior storage
    ! -- Classes
    type, abstract :: sparse_t
      integer :: base = 1 !! index base = 0 for (C) or 1 (Fortran)
      integer :: nrows= 0 !! number of rows
      integer :: ncols= 0 !! number of columns
      integer :: NNZ  = 0 !! number of non-zero values
      integer :: sym  = k_NOSYMMETRY !! assumed storage symmetry
    end type
  
    !! COO: COOrdinates compresed format
    type, extends(sparse_t) :: COO_t
      logical               :: isOrdered = .false. !! wether the matrix is ordered or not
      integer, allocatable  :: index(:,:) !! Matrix coordinates index(2,nnz)
    contains
      procedure :: malloc => malloc_coo
    end type
  
    type, extends(COO_t) :: COOr32_t
      real(sp), allocatable :: data(:) !! single precision values
    end type
  
    type, extends(COO_t) :: COOr64_t
      real(dp), allocatable :: data(:) !! double precision values
    end type
   
    type, extends(COO_t) :: COOc32_t
      complex(sp), allocatable :: data(:) !! single precision values
    end type
  
    type, extends(COO_t) :: COOc64_t
      complex(dp), allocatable :: data(:) !! double precision values
    end type
  
    !! CSR: Compressed sparse row or Yale format
    type, extends(sparse_t) :: CSR_t  
      integer, allocatable  :: col(:)    !! matrix column pointer
      integer, allocatable  :: rowptr(:) !! matrix row pointer
    contains
      procedure :: malloc => malloc_csr
    end type
  
    type, extends(CSR_t) :: CSRr32_t
      real(sp), allocatable :: data(:)   !! single precision values
    end type
  
    type, extends(CSR_t) :: CSRr64_t 
      real(dp), allocatable :: data(:)   !! double precision values
    end type
  
    type, extends(CSR_t) :: CSRc32_t
      complex(sp), allocatable :: data(:)   !! single precision values
    end type
  
    type, extends(CSR_t) :: CSRc64_t 
      complex(dp), allocatable :: data(:)   !! double precision values
    end type
  

    !! CSC: Compressed sparse column
    type, extends(sparse_t) :: CSC_t  
      integer, allocatable  :: colptr(:) !! matrix column pointer
      integer, allocatable  :: row(:)    !! matrix row pointer
    contains
      procedure :: malloc => malloc_csc
    end type
  
    type, extends(CSC_t) :: CSCr32_t 
      real(sp), allocatable :: data(:)   !! single precision values
    end type
  
    type, extends(CSC_t) :: CSCr64_t 
      real(dp), allocatable :: data(:)   !! double precision values
    end type
  
    type, extends(CSC_t) :: CSCc32_t 
      complex(sp), allocatable :: data(:)   !! single precision values
    end type
  
    type, extends(CSC_t) :: CSCc64_t 
      complex(dp), allocatable :: data(:)   !! double precision values
    end type
  
    !! Compressed ELLPACK
    type, extends(sparse_t) :: ELL_t 
      integer               :: K = 0 !! maximum number of nonzeros per row
      integer, allocatable  :: index(:,:) !! column indices
    contains
      procedure :: malloc => malloc_ell
    end type
  
    type, extends(ELL_t) :: ELLr32_t 
      real(sp), allocatable :: data(:,:) !! single precision values
    end type
  
    type, extends(ELL_t) :: ELLr64_t 
      real(dp), allocatable :: data(:,:) !! double precision values
    end type
  
  
    type, extends(ELL_t) :: ELLc32_t 
      complex(sp), allocatable :: data(:,:) !! single precision values
    end type
  
    type, extends(ELL_t) :: ELLc64_t 
      complex(dp), allocatable :: data(:,:) !! double precision values
    end type
  
  contains
  
  subroutine malloc_coo(self,num_rows,num_cols,nnz)
      class(COO_t) :: self
      integer, intent(in) :: num_rows
      integer, intent(in) :: num_cols
      integer, intent(in) :: nnz
  
      integer,  allocatable :: itemp(:,:)
      real(sp), allocatable :: stemp(:)
      real(dp), allocatable :: dtemp(:)
      complex(sp),allocatable:: ctemp(:)
      complex(dp),allocatable:: ztemp(:)
      
      self%nrows = num_rows
      self%ncols = num_cols
      self%NNZ = nnz
  
      if(.not.allocated(self%index)) then
          allocate(itemp(2,nnz) , source = 0 )
      else
          allocate(itemp(2,nnz) , source = self%index )
      end if
      call move_alloc(from=itemp,to=self%index)
  
      select type(self)
          type is(COOr32_t)
              if(.not.allocated(self%data)) then
                  allocate(stemp(nnz) , source = 0._sp )
              else
                  allocate(stemp(nnz) , source = self%data )
              end if
              call move_alloc(from=stemp,to=self%data)
          type is(COOr64_t)
              if(.not.allocated(self%data)) then
                  allocate(dtemp(nnz) , source = 0._dp )
              else
                  allocate(dtemp(nnz) , source = self%data )
              end if
              call move_alloc(from=dtemp,to=self%data)
          type is(COOc32_t)
              if(.not.allocated(self%data)) then
                  allocate(ctemp(nnz) , source = (0._sp,0._sp) )
              else
                  allocate(ctemp(nnz) , source = self%data )
              end if
              call move_alloc(from=ctemp,to=self%data)
          type is(COOc64_t)
              if(.not.allocated(self%data)) then
                  allocate(ztemp(nnz) , source = (0._dp,0._dp) )
              else
                  allocate(ztemp(nnz) , source = self%data )
              end if
              call move_alloc(from=ztemp,to=self%data)
      end select
  end subroutine
  
  subroutine malloc_csr(self,num_rows,num_cols,nnz)
      class(CSR_t) :: self
      integer, intent(in) :: num_rows
      integer, intent(in) :: num_cols
      integer, intent(in) :: nnz
  
      integer,  allocatable :: itemp(:)
      real(sp), allocatable :: stemp(:)
      real(dp), allocatable :: dtemp(:)
      complex(sp),allocatable:: ctemp(:)
      complex(dp),allocatable:: ztemp(:)
  
      self%nrows = num_rows
      self%ncols = num_cols
      self%NNZ = nnz
  
      if(.not.allocated(self%col)) then
          allocate(itemp(nnz) , source = 0 )
      else
          allocate(itemp(nnz) , source = self%col )
      end if
      call move_alloc(from=itemp,to=self%col)
  
      if(.not.allocated(self%rowptr)) then
          allocate(itemp(num_rows+1) , source = 0 )
      else
          allocate(itemp(num_rows+1) , source = self%rowptr )
      end if
      call move_alloc(from=itemp,to=self%rowptr)
  
      select type(self)
          type is(CSRr32_t)
              if(.not.allocated(self%data)) then
                  allocate(stemp(nnz) , source = 0._sp )
              else
                  allocate(stemp(nnz) , source = self%data )
              end if
              call move_alloc(from=stemp,to=self%data)
          type is(CSRr64_t)
              if(.not.allocated(self%data)) then
                  allocate(dtemp(nnz) , source = 0._dp )
              else
                  allocate(dtemp(nnz) , source = self%data )
              end if
              call move_alloc(from=dtemp,to=self%data)
          type is(CSRc32_t)
              if(.not.allocated(self%data)) then
                  allocate(ctemp(nnz) , source = (0._sp,0._sp) )
              else
                  allocate(ctemp(nnz) , source = self%data )
              end if
              call move_alloc(from=ctemp,to=self%data)
          type is(CSRc64_t)
              if(.not.allocated(self%data)) then
                  allocate(ztemp(nnz) , source = (0._dp,0._dp) )
              else
                  allocate(ztemp(nnz) , source = self%data )
              end if
              call move_alloc(from=ztemp,to=self%data)
      end select
  end subroutine
  
  subroutine malloc_csc(self,num_rows,num_cols,nnz)
      class(CSC_t) :: self
      integer, intent(in) :: num_rows
      integer, intent(in) :: num_cols
      integer, intent(in) :: nnz
  
      integer,  allocatable :: itemp(:)
      real(sp), allocatable :: stemp(:)
      real(dp), allocatable :: dtemp(:)
      complex(sp),allocatable:: ctemp(:)
      complex(dp),allocatable:: ztemp(:)
  
      self%nrows = num_rows
      self%ncols = num_cols
      self%NNZ = nnz
  
      if(.not.allocated(self%row)) then
          allocate(itemp(nnz) , source = 0 )
      else
          allocate(itemp(nnz) , source = self%row )
      end if
      call move_alloc(from=itemp,to=self%row)
  
      if(.not.allocated(self%colptr)) then
          allocate(itemp(num_cols+1) , source = 0 )
      else
          allocate(itemp(num_cols+1) , source = self%colptr )
      end if
      call move_alloc(from=itemp,to=self%colptr)
  
      select type(self)
          type is(CSCr32_t)
              if(.not.allocated(self%data)) then
                  allocate(stemp(nnz) , source = 0._sp )
              else
                  allocate(stemp(nnz) , source = self%data )
              end if
              call move_alloc(from=stemp,to=self%data)
          type is(CSCr64_t)
              if(.not.allocated(self%data)) then
                  allocate(dtemp(nnz) , source = 0._dp )
              else
                  allocate(dtemp(nnz) , source = self%data )
              end if
              call move_alloc(from=dtemp,to=self%data)
          type is(CSCc32_t)
              if(.not.allocated(self%data)) then
                  allocate(ctemp(nnz) , source = (0._sp,0._sp) )
              else
                  allocate(ctemp(nnz) , source = self%data )
              end if
              call move_alloc(from=ctemp,to=self%data)
          type is(CSCc64_t)
              if(.not.allocated(self%data)) then
                  allocate(ztemp(nnz) , source = (0._dp,0._dp) )
              else
                  allocate(ztemp(nnz) , source = self%data )
              end if
              call move_alloc(from=ztemp,to=self%data)
      end select
  end subroutine
  
  subroutine malloc_ell(self,num_rows,num_cols,num_nz_rows)
      class(ELL_t) :: self
      integer, intent(in) :: num_rows    !! number of rows
      integer, intent(in) :: num_cols    !! number of columns
      integer, intent(in) :: num_nz_rows !! number of non zeros per row
  
      integer,  allocatable :: itemp(:,:)
      real(sp), allocatable :: stemp(:,:)
      real(dp), allocatable :: dtemp(:,:)
      complex(sp),allocatable:: ctemp(:,:)
      complex(dp),allocatable:: ztemp(:,:)
  
      self%nrows = num_rows
      self%ncols = num_cols
      self%K = num_nz_rows
  
      if(.not.allocated(self%index)) then
          allocate(itemp(num_rows,num_nz_rows) , source = 0 )
      else
          allocate(itemp(num_rows,num_nz_rows) , source = self%index )
      end if
      call move_alloc(from=itemp,to=self%index)
  
      select type(self)
          type is(ELLr32_t)
              if(.not.allocated(self%data)) then
                  allocate(stemp(num_rows,num_nz_rows) , source = 0._sp )
              else
                  allocate(stemp(num_rows,num_nz_rows) , source = self%data )
              end if
              call move_alloc(from=stemp,to=self%data)
          type is(ELLr64_t)
              if(.not.allocated(self%data)) then
                  allocate(dtemp(num_rows,num_nz_rows) , source = 0._dp )
              else
                  allocate(dtemp(num_rows,num_nz_rows) , source = self%data )
              end if
              call move_alloc(from=dtemp,to=self%data)
          type is(ELLc32_t)
              if(.not.allocated(self%data)) then
                  allocate(ctemp(num_rows,num_nz_rows) , source = (0._sp,0._sp) )
              else
                  allocate(ctemp(num_rows,num_nz_rows) , source = self%data )
              end if
              call move_alloc(from=ctemp,to=self%data)
          type is(ELLc64_t)
              if(.not.allocated(self%data)) then
                  allocate(ztemp(num_rows,num_nz_rows) , source = (0._dp,0._dp) )
              else
                  allocate(ztemp(num_rows,num_nz_rows) , source = self%data )
              end if
              call move_alloc(from=ztemp,to=self%data)
      end select
  end subroutine

  end module fsparse_matrix_gallery
