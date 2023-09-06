!---------------------------------------------------
! Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module conversions
    use iso_fortran_env
    use matrix_gallery
    use sparse_sort
    implicit none
    
    interface dense2coo
        module procedure dense2coo_sp
        module procedure dense2coo_dp
    end interface

    interface coo2dense
        module procedure coo2dense_sp
        module procedure coo2dense_dp
    end interface

    interface coo2csr
        module procedure coo2csr_ordered
        module procedure coo2csr_ordered_sp
        module procedure coo2csr_ordered_dp
    end interface

contains

    subroutine dense2coo_sp(dense,COO)
        integer, parameter :: wp = real32
        real(wp), intent(in) :: dense(:,:)
        type(COOr32_t), intent(inout) :: COO
        integer :: num_rows, num_cols, nnz
        integer :: i, j, idx

        num_rows = size(dense,dim=1)
        num_cols = size(dense,dim=2)
        nnz = count( abs(dense) > tiny(1._wp) )

        call COO%malloc(num_rows,num_cols,nnz)

        idx = 1
        do i = 1, num_rows
            do j = 1, num_cols
                if(abs(dense(i,j)) < tiny(1._wp)) cycle
                COO%index(1,idx) = i
                COO%index(2,idx) = j
                COO%data(idx) = dense(i,j)
                idx = idx + 1
            end do
        end do
        COO%isOrdered = .true.
    end subroutine

    subroutine dense2coo_dp(dense,COO)
        integer, parameter :: wp = real64
        real(wp), intent(in) :: dense(:,:)
        type(COOr64_t), intent(inout) :: COO
        integer :: num_rows, num_cols, nnz
        integer :: i, j, idx

        num_rows = size(dense,dim=1)
        num_cols = size(dense,dim=2)
        nnz = count( abs(dense) > tiny(1._wp) )

        call COO%malloc(num_rows,num_cols,nnz)

        idx = 1
        do i = 1, num_rows
            do j = 1, num_cols
                if(abs(dense(i,j)) < tiny(1._wp)) cycle
                COO%index(1,idx) = i
                COO%index(2,idx) = j
                COO%data(idx) = dense(i,j)
                idx = idx + 1
            end do
        end do
        COO%isOrdered = .true.
    end subroutine

    subroutine coo2dense_sp(COO,dense)
        integer, parameter :: wp = real32
        type(COOr32_t), intent(in) :: COO
        real(wp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2dense_dp(COO,dense)
        integer, parameter :: wp = real64
        type(COOr64_t), intent(in) :: COO
        real(wp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2csr_ordered(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COO_t), intent(in)    :: COO
        type(CSR_t), intent(inout) :: CSR
        integer :: i

        associate( nnz=>COO%nnz, num_rows=>COO%N, num_cols=>COO%M, base=>COO%base )
        CSR%NNZ = nnz; CSR%N = num_rows; CSR%M = num_cols; CSR%base = base

        if( allocated(CSR%col) ) then
            CSR%col(1:nnz)  = COO%index(2,1:nnz)
            CSR%rowptr(1:num_rows) = 0
        else 
            allocate( CSR%col(nnz)  , source = COO%index(2,1:nnz) )
            allocate( CSR%rowptr(num_rows+1) , source = 0 )
        end if

        CSR%rowptr(1) = 1
        do i = 1, nnz
            CSR%rowptr( COO%index(1,i)+1 ) = CSR%rowptr( COO%index(1,i)+1 ) + 1
        end do
        do i = 1, num_rows
            CSR%rowptr( i+1 ) = CSR%rowptr( i+1 ) + CSR%rowptr( i )
        end do
        end associate
    end subroutine

    subroutine coo2csr_ordered_sp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COOr32_t), intent(in)    :: COO
        type(CSRr32_t), intent(inout) :: CSR
        integer :: i

        associate( nnz=>COO%nnz, num_rows=>COO%N, num_cols=>COO%M, base=>COO%base )
        CSR%NNZ = nnz; CSR%N = num_rows; CSR%M = num_cols; CSR%base = base
        if( allocated(CSR%col) ) then
            CSR%col(1:nnz)  = COO%index(2,1:nnz)
            CSR%rowptr(1:num_rows) = 0
            CSR%data(1:nnz) = COO%data(1:nnz)
        else 
            allocate( CSR%col(nnz)  , source = COO%index(2,1:nnz) )
            allocate( CSR%rowptr(num_rows+1) , source = 0 )
            allocate( CSR%data(nnz) , source = COO%data(1:nnz) )
        end if

        CSR%rowptr(1) = 1
        do i = 1, nnz
            CSR%rowptr( COO%index(1,i)+1 ) = CSR%rowptr( COO%index(1,i)+1 ) + 1
        end do
        do i = 1, num_rows
            CSR%rowptr( i+1 ) = CSR%rowptr( i+1 ) + CSR%rowptr( i )
        end do
        end associate
    end subroutine

    subroutine coo2csr_ordered_dp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COOr64_t), intent(in)    :: COO
        type(CSRr64_t), intent(inout) :: CSR
        integer :: i

        associate( nnz=>COO%nnz, num_rows=>COO%N, num_cols=>COO%M, base=>COO%base )
        CSR%NNZ = nnz; CSR%N = num_rows; CSR%M = num_cols; CSR%base = base
        if( allocated(CSR%col) ) then
            CSR%col(1:nnz)  = COO%index(2,1:nnz)
            CSR%rowptr(1:num_rows) = 0
            CSR%data(1:nnz) = COO%data(1:nnz)
        else 
            allocate( CSR%col(nnz)  , source = COO%index(2,1:nnz) )
            allocate( CSR%rowptr(num_rows+1) , source = 0 )
            allocate( CSR%data(nnz) , source = COO%data(1:nnz) )
        end if

        CSR%rowptr(1) = 1
        do i = 1, nnz
            CSR%rowptr( COO%index(1,i)+1 ) = CSR%rowptr( COO%index(1,i)+1 ) + 1
        end do
        do i = 1, num_rows
            CSR%rowptr( i+1 ) = CSR%rowptr( i+1 ) + CSR%rowptr( i )
        end do
        end associate
    end subroutine
    
end module conversions