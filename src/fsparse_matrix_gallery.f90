!>---------------------------------------------------
!> Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!>
!> Use of this source code is governed by a MIT
!> license that can be found in the LICENSE.md file
!>---------------------------------------------------
module fsparse_matrix_gallery
    use fsparse_constants
    implicit none
    
    private
  
    ! -- Global parameters
    enum, bind(C)
        enumerator :: k_NOSYMMETRY !! Full Sparse matrix (no symmetry considerations)
        enumerator :: k_SYMTRIINF  !! Symmetric Sparse matrix with triangular inferior storage
        enumerator :: k_SYMTRISUP  !! Symmetric Sparse matrix with triangular supperior storage
    end enum
    public :: k_NOSYMMETRY, k_SYMTRIINF, k_SYMTRISUP
    ! -- Classes
    type, public, abstract :: sparse_t
      integer :: nrows = 0 !! number of rows
      integer :: ncols = 0 !! number of columns
      integer :: nnz   = 0 !! number of non-zero values
      integer :: sym   = k_NOSYMMETRY !! assumed storage symmetry
      integer :: base  = 1 !! index base = 0 for (C) or 1 (Fortran)
    end type

    !! COO: COOrdinates compresed format
    type, public, extends(sparse_t) :: COO_t
      logical               :: isOrdered = .false. !! wether the matrix is ordered or not
      integer, allocatable  :: index(:,:) !! Matrix coordinates index(2,nnz)
    contains
      procedure :: malloc => malloc_coo
    end type

    type, public, extends(COO_t) :: COO_sp
        real(sp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_coo_sp
        procedure :: set => set_value_coo_sp
    end type
    type, public, extends(COO_t) :: COO_dp
        real(dp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_coo_dp
        procedure :: set => set_value_coo_dp
    end type
    type, public, extends(COO_t) :: COO_qp
        real(qp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_coo_qp
        procedure :: set => set_value_coo_qp
    end type
    type, public, extends(COO_t) :: COO_csp
        complex(sp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_coo_csp
        procedure :: set => set_value_coo_csp
    end type
    type, public, extends(COO_t) :: COO_cdp
        complex(dp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_coo_cdp
        procedure :: set => set_value_coo_cdp
    end type
    type, public, extends(COO_t) :: COO_cqp
        complex(qp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_coo_cqp
        procedure :: set => set_value_coo_cqp
    end type

    !! CSR: Compressed sparse row or Yale format
    type, public, extends(sparse_t) :: CSR_t  
      integer, allocatable  :: col(:)    !! matrix column pointer
      integer, allocatable  :: rowptr(:) !! matrix row pointer
    contains
      procedure :: malloc => malloc_csr
    end type
  
    type, public, extends(CSR_t) :: CSR_sp
        real(sp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csr_sp
        procedure :: set => set_value_csr_sp
    end type
    type, public, extends(CSR_t) :: CSR_dp
        real(dp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csr_dp
        procedure :: set => set_value_csr_dp
    end type
    type, public, extends(CSR_t) :: CSR_qp
        real(qp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csr_qp
        procedure :: set => set_value_csr_qp
    end type
    type, public, extends(CSR_t) :: CSR_csp
        complex(sp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csr_csp
        procedure :: set => set_value_csr_csp
    end type
    type, public, extends(CSR_t) :: CSR_cdp
        complex(dp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csr_cdp
        procedure :: set => set_value_csr_cdp
    end type
    type, public, extends(CSR_t) :: CSR_cqp
        complex(qp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csr_cqp
        procedure :: set => set_value_csr_cqp
    end type

    !! CSC: Compressed sparse column
    type, public, extends(sparse_t) :: CSC_t  
      integer, allocatable  :: colptr(:) !! matrix column pointer
      integer, allocatable  :: row(:)    !! matrix row pointer
    contains
      procedure :: malloc => malloc_csc
    end type
  
    type, public, extends(CSC_t) :: CSC_sp
        real(sp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csc_sp
        procedure :: set => set_value_csc_sp
    end type
    type, public, extends(CSC_t) :: CSC_dp
        real(dp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csc_dp
        procedure :: set => set_value_csc_dp
    end type
    type, public, extends(CSC_t) :: CSC_qp
        real(qp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csc_qp
        procedure :: set => set_value_csc_qp
    end type
    type, public, extends(CSC_t) :: CSC_csp
        complex(sp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csc_csp
        procedure :: set => set_value_csc_csp
    end type
    type, public, extends(CSC_t) :: CSC_cdp
        complex(dp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csc_cdp
        procedure :: set => set_value_csc_cdp
    end type
    type, public, extends(CSC_t) :: CSC_cqp
        complex(qp), allocatable :: data(:) 
    contains
        procedure :: get => get_value_csc_cqp
        procedure :: set => set_value_csc_cqp
    end type
  
    !! Compressed ELLPACK
    type, public, extends(sparse_t) :: ELL_t 
      integer               :: K = 0 !! maximum number of nonzeros per row
      integer, allocatable  :: index(:,:) !! column indices
    contains
      procedure :: malloc => malloc_ell
    end type
  
    type, public, extends(ELL_t) :: ELL_sp
        real(sp), allocatable :: data(:,:) 
    contains
        procedure :: get => get_value_ell_sp
        procedure :: set => set_value_ell_sp
    end type
    type, public, extends(ELL_t) :: ELL_dp
        real(dp), allocatable :: data(:,:) 
    contains
        procedure :: get => get_value_ell_dp
        procedure :: set => set_value_ell_dp
    end type
    type, public, extends(ELL_t) :: ELL_qp
        real(qp), allocatable :: data(:,:) 
    contains
        procedure :: get => get_value_ell_qp
        procedure :: set => set_value_ell_qp
    end type
    type, public, extends(ELL_t) :: ELL_csp
        complex(sp), allocatable :: data(:,:) 
    contains
        procedure :: get => get_value_ell_csp
        procedure :: set => set_value_ell_csp
    end type
    type, public, extends(ELL_t) :: ELL_cdp
        complex(dp), allocatable :: data(:,:) 
    contains
        procedure :: get => get_value_ell_cdp
        procedure :: set => set_value_ell_cdp
    end type
    type, public, extends(ELL_t) :: ELL_cqp
        complex(qp), allocatable :: data(:,:) 
    contains
        procedure :: get => get_value_ell_cqp
        procedure :: set => set_value_ell_cqp
    end type

contains

    subroutine malloc_coo(self,num_rows,num_cols,nnz)
        class(COO_t) :: self
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols
        integer, intent(in) :: nnz

        integer,  allocatable :: temp_idx(:,:)
        !-----------------------------------------------------

        self%nrows = num_rows
        self%ncols = num_cols
        self%nnz   = nnz

        if(.not.allocated(self%index)) then
            allocate(temp_idx(2,nnz) , source = 0 )
        else
            allocate(temp_idx(2,nnz) , source = self%index )
        end if
        call move_alloc(from=temp_idx,to=self%index)

        select type(self)
            type is(COO_sp)
                block
                real(sp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_sp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(COO_dp)
                block
                real(dp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_dp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(COO_qp)
                block
                real(qp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_qp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(COO_csp)
                block
                complex(sp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_csp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(COO_cdp)
                block
                complex(dp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_cdp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(COO_cqp)
                block
                complex(qp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_cqp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
        end select
    end subroutine

    subroutine malloc_csr(self,num_rows,num_cols,nnz)
        class(CSR_t) :: self
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols
        integer, intent(in) :: nnz

        integer,  allocatable :: temp_idx(:)
        !-----------------------------------------------------

        self%nrows = num_rows
        self%ncols = num_cols
        self%nnz   = nnz

        if(.not.allocated(self%col)) then
            allocate(temp_idx(nnz) , source = 0 )
        else
            allocate(temp_idx(nnz) , source = self%col )
        end if
        call move_alloc(from=temp_idx,to=self%col)

        if(.not.allocated(self%rowptr)) then
            allocate(temp_idx(num_rows+1) , source = 0 )
        else
            allocate(temp_idx(num_rows+1) , source = self%rowptr )
        end if
        call move_alloc(from=temp_idx,to=self%rowptr)

        select type(self)
            type is(CSR_sp)
                block
                real(sp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_sp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSR_dp)
                block
                real(dp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_dp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSR_qp)
                block
                real(qp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_qp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSR_csp)
                block
                complex(sp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_csp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSR_cdp)
                block
                complex(dp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_cdp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSR_cqp)
                block
                complex(qp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_cqp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
        end select
    end subroutine

    subroutine malloc_csc(self,num_rows,num_cols,nnz)
        class(CSC_t) :: self
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols
        integer, intent(in) :: nnz

        integer,  allocatable :: temp_idx(:)
        !-----------------------------------------------------

        self%nrows = num_rows
        self%ncols = num_cols
        self%nnz   = nnz

        if(.not.allocated(self%row)) then
            allocate(temp_idx(nnz) , source = 0 )
        else
            allocate(temp_idx(nnz) , source = self%row )
        end if
        call move_alloc(from=temp_idx,to=self%row)

        if(.not.allocated(self%colptr)) then
            allocate(temp_idx(num_cols+1) , source = 0 )
        else
            allocate(temp_idx(num_cols+1) , source = self%colptr )
        end if
        call move_alloc(from=temp_idx,to=self%colptr)

        select type(self)
            type is(CSC_sp)
                block
                real(sp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_sp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSC_dp)
                block
                real(dp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_dp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSC_qp)
                block
                real(qp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_qp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSC_csp)
                block
                complex(sp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_csp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSC_cdp)
                block
                complex(dp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_cdp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(CSC_cqp)
                block
                complex(qp), allocatable :: temp(:)
                if(.not.allocated(self%data)) then
                    allocate(temp(nnz)); temp = zero_cqp
                else
                    allocate(temp(nnz) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
        end select
    end subroutine

    subroutine malloc_ell(self,num_rows,num_cols,num_nz_rows)
        class(ELL_t) :: self
        integer, intent(in) :: num_rows    !! number of rows
        integer, intent(in) :: num_cols    !! number of columns
        integer, intent(in) :: num_nz_rows !! number of non zeros per row

        integer,  allocatable :: temp_idx(:,:)
        !-----------------------------------------------------

        self%nrows = num_rows
        self%ncols = num_cols
        self%K     = num_nz_rows

        if(.not.allocated(self%index)) then
            allocate(temp_idx(num_rows,num_nz_rows) , source = 0 )
        else
            allocate(temp_idx(num_rows,num_nz_rows) , source = self%index )
        end if
        call move_alloc(from=temp_idx,to=self%index)

        select type(self)
            type is(ELL_sp)
                block
                real(sp), allocatable :: temp(:,:)
                if(.not.allocated(self%data)) then
                    allocate(temp(num_rows,num_nz_rows)); temp = zero_sp
                else
                    allocate(temp(num_rows,num_nz_rows) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(ELL_dp)
                block
                real(dp), allocatable :: temp(:,:)
                if(.not.allocated(self%data)) then
                    allocate(temp(num_rows,num_nz_rows)); temp = zero_dp
                else
                    allocate(temp(num_rows,num_nz_rows) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(ELL_qp)
                block
                real(qp), allocatable :: temp(:,:)
                if(.not.allocated(self%data)) then
                    allocate(temp(num_rows,num_nz_rows)); temp = zero_qp
                else
                    allocate(temp(num_rows,num_nz_rows) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(ELL_csp)
                block
                complex(sp), allocatable :: temp(:,:)
                if(.not.allocated(self%data)) then
                    allocate(temp(num_rows,num_nz_rows)); temp = zero_csp
                else
                    allocate(temp(num_rows,num_nz_rows) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(ELL_cdp)
                block
                complex(dp), allocatable :: temp(:,:)
                if(.not.allocated(self%data)) then
                    allocate(temp(num_rows,num_nz_rows)); temp = zero_cdp
                else
                    allocate(temp(num_rows,num_nz_rows) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
            type is(ELL_cqp)
                block
                complex(qp), allocatable :: temp(:,:)
                if(.not.allocated(self%data)) then
                    allocate(temp(num_rows,num_nz_rows)); temp = zero_cqp
                else
                    allocate(temp(num_rows,num_nz_rows) , source = self%data )
                end if
                call move_alloc(from=temp,to=self%data)
                end block
        end select
    end subroutine

    !==================================================================
    ! data accessors
    !==================================================================

    pure subroutine get_value_coo_sp(self,val,ik,jk)
        class(COO_sp), intent(in) :: self
        real(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_sp
    end subroutine

    subroutine set_value_coo_sp(self,val,ik,jk)
        class(COO_sp), intent(inout) :: self
        real(sp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_coo_dp(self,val,ik,jk)
        class(COO_dp), intent(in) :: self
        real(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_dp
    end subroutine

    subroutine set_value_coo_dp(self,val,ik,jk)
        class(COO_dp), intent(inout) :: self
        real(dp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_coo_qp(self,val,ik,jk)
        class(COO_qp), intent(in) :: self
        real(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_qp
    end subroutine

    subroutine set_value_coo_qp(self,val,ik,jk)
        class(COO_qp), intent(inout) :: self
        real(qp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_coo_csp(self,val,ik,jk)
        class(COO_csp), intent(in) :: self
        complex(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_csp
    end subroutine

    subroutine set_value_coo_csp(self,val,ik,jk)
        class(COO_csp), intent(inout) :: self
        complex(sp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_coo_cdp(self,val,ik,jk)
        class(COO_cdp), intent(in) :: self
        complex(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_cdp
    end subroutine

    subroutine set_value_coo_cdp(self,val,ik,jk)
        class(COO_cdp), intent(inout) :: self
        complex(dp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_coo_cqp(self,val,ik,jk)
        class(COO_cqp), intent(in) :: self
        complex(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_cqp
    end subroutine

    subroutine set_value_coo_cqp(self,val,ik,jk)
        class(COO_cqp), intent(inout) :: self
        complex(qp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1,self%nnz
            if( ik == self%index(1,k) .and. jk == self%index(2,k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine


    pure subroutine get_value_csr_sp(self,val,ik,jk)
        class(CSR_sp), intent(in) :: self
        real(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_sp
    end subroutine

    subroutine set_value_csr_sp(self,val,ik,jk)
        class(CSR_sp), intent(inout) :: self
        real(sp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csr_dp(self,val,ik,jk)
        class(CSR_dp), intent(in) :: self
        real(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_dp
    end subroutine

    subroutine set_value_csr_dp(self,val,ik,jk)
        class(CSR_dp), intent(inout) :: self
        real(dp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csr_qp(self,val,ik,jk)
        class(CSR_qp), intent(in) :: self
        real(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_qp
    end subroutine

    subroutine set_value_csr_qp(self,val,ik,jk)
        class(CSR_qp), intent(inout) :: self
        real(qp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csr_csp(self,val,ik,jk)
        class(CSR_csp), intent(in) :: self
        complex(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_csp
    end subroutine

    subroutine set_value_csr_csp(self,val,ik,jk)
        class(CSR_csp), intent(inout) :: self
        complex(sp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csr_cdp(self,val,ik,jk)
        class(CSR_cdp), intent(in) :: self
        complex(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_cdp
    end subroutine

    subroutine set_value_csr_cdp(self,val,ik,jk)
        class(CSR_cdp), intent(inout) :: self
        complex(dp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csr_cqp(self,val,ik,jk)
        class(CSR_cqp), intent(in) :: self
        complex(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_cqp
    end subroutine

    subroutine set_value_csr_cqp(self,val,ik,jk)
        class(CSR_cqp), intent(inout) :: self
        complex(qp), intent(in) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%rowptr(ik), self%rowptr(ik+1)-1
            if( jk == self%col(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine


    pure subroutine get_value_csc_sp(self,val,ik,jk)
        class(CSC_sp), intent(in) :: self
        real(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_sp
    end subroutine

    subroutine set_value_csc_sp(self,val,ik,jk)
        class(CSC_sp), intent(inout) :: self
        real(sp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csc_dp(self,val,ik,jk)
        class(CSC_dp), intent(in) :: self
        real(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_dp
    end subroutine

    subroutine set_value_csc_dp(self,val,ik,jk)
        class(CSC_dp), intent(inout) :: self
        real(dp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csc_qp(self,val,ik,jk)
        class(CSC_qp), intent(in) :: self
        real(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_qp
    end subroutine

    subroutine set_value_csc_qp(self,val,ik,jk)
        class(CSC_qp), intent(inout) :: self
        real(qp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csc_csp(self,val,ik,jk)
        class(CSC_csp), intent(in) :: self
        complex(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_csp
    end subroutine

    subroutine set_value_csc_csp(self,val,ik,jk)
        class(CSC_csp), intent(inout) :: self
        complex(sp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csc_cdp(self,val,ik,jk)
        class(CSC_cdp), intent(in) :: self
        complex(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_cdp
    end subroutine

    subroutine set_value_csc_cdp(self,val,ik,jk)
        class(CSC_cdp), intent(inout) :: self
        complex(dp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_csc_cqp(self,val,ik,jk)
        class(CSC_cqp), intent(in) :: self
        complex(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                val = self%data(k)
                return
            end if
        end do
        val = zero_cqp
    end subroutine

    subroutine set_value_csc_cqp(self,val,ik,jk)
        class(CSC_cqp), intent(inout) :: self
        complex(qp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = self%colptr(jk), self%colptr(jk+1)-1
            if( ik == self%row(k) ) then
                self%data(k) = val
                return
            end if
        end do
    end subroutine


    pure subroutine get_value_ell_sp(self,val,ik,jk)
        class(ELL_sp), intent(in) :: self
        real(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                val = self%data(ik,k)
                return
            end if
        end do
        val = zero_sp
    end subroutine

    subroutine set_value_ell_sp(self,val,ik,jk)
        class(ELL_sp), intent(inout) :: self
        real(sp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                self%data(ik,k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_ell_dp(self,val,ik,jk)
        class(ELL_dp), intent(in) :: self
        real(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                val = self%data(ik,k)
                return
            end if
        end do
        val = zero_dp
    end subroutine

    subroutine set_value_ell_dp(self,val,ik,jk)
        class(ELL_dp), intent(inout) :: self
        real(dp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                self%data(ik,k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_ell_qp(self,val,ik,jk)
        class(ELL_qp), intent(in) :: self
        real(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                val = self%data(ik,k)
                return
            end if
        end do
        val = zero_qp
    end subroutine

    subroutine set_value_ell_qp(self,val,ik,jk)
        class(ELL_qp), intent(inout) :: self
        real(qp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                self%data(ik,k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_ell_csp(self,val,ik,jk)
        class(ELL_csp), intent(in) :: self
        complex(sp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                val = self%data(ik,k)
                return
            end if
        end do
        val = zero_csp
    end subroutine

    subroutine set_value_ell_csp(self,val,ik,jk)
        class(ELL_csp), intent(inout) :: self
        complex(sp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                self%data(ik,k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_ell_cdp(self,val,ik,jk)
        class(ELL_cdp), intent(in) :: self
        complex(dp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                val = self%data(ik,k)
                return
            end if
        end do
        val = zero_cdp
    end subroutine

    subroutine set_value_ell_cdp(self,val,ik,jk)
        class(ELL_cdp), intent(inout) :: self
        complex(dp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                self%data(ik,k) = val
                return
            end if
        end do
    end subroutine

    pure subroutine get_value_ell_cqp(self,val,ik,jk)
        class(ELL_cqp), intent(in) :: self
        complex(qp), intent(out) :: val
        integer, intent(in) :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                val = self%data(ik,k)
                return
            end if
        end do
        val = zero_cqp
    end subroutine

    subroutine set_value_ell_cqp(self,val,ik,jk)
        class(ELL_cqp), intent(inout) :: self
        complex(qp), intent(in) :: val
        integer, intent(in)  :: ik, jk
        integer :: k
        ! naive implementation
        do k = 1 , self%K
            if( jk == self%index(ik,k) ) then
                self%data(ik,k) = val
                return
            end if
        end do
    end subroutine

    
end module fsparse_matrix_gallery