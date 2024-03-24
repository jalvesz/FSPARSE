!---------------------------------------------------
! Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module fsparse_conversions
    use iso_fortran_env
    use fsparse_matrix_gallery
    use fsparse_sort
    implicit none
    
    interface dense2coo
        module procedure dense2coo_r_sp
        module procedure dense2coo_r_dp
        module procedure dense2coo_c_sp
        module procedure dense2coo_c_dp
    end interface

    interface coo2dense
        module procedure coo2dense_r_sp
        module procedure coo2dense_r_dp
        module procedure coo2dense_c_sp
        module procedure coo2dense_c_dp
    end interface

    interface coo2csr
        module procedure coo2csr_ordered
        module procedure coo2csr_ordered_r_sp
        module procedure coo2csr_ordered_r_dp
        module procedure coo2csr_ordered_c_sp
        module procedure coo2csr_ordered_c_dp
    end interface

    interface csr2coo
        module procedure csr2coo
        module procedure csr2coo_r_sp
        module procedure csr2coo_r_dp
        module procedure csr2coo_c_sp
        module procedure csr2coo_c_dp
    end interface

contains

    subroutine dense2coo_r_sp(dense,COO)
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

    subroutine dense2coo_r_dp(dense,COO)
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

    subroutine dense2coo_c_sp(dense,COO)
        integer, parameter :: wp = real32
        complex(wp), intent(in) :: dense(:,:)
        type(COOc32_t), intent(inout) :: COO
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

    subroutine dense2coo_c_dp(dense,COO)
        integer, parameter :: wp = real64
        complex(wp), intent(in) :: dense(:,:)
        type(COOc64_t), intent(inout) :: COO
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
    
    subroutine coo2dense_r_sp(COO,dense)
        integer, parameter :: wp = real32
        type(COOr32_t), intent(in) :: COO
        real(wp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2dense_r_dp(COO,dense)
        integer, parameter :: wp = real64
        type(COOr64_t), intent(in) :: COO
        real(wp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine
    
    subroutine coo2dense_c_sp(COO,dense)
        integer, parameter :: wp = real32
        type(COOc32_t), intent(in) :: COO
        complex(wp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2dense_c_dp(COO,dense)
        integer, parameter :: wp = real64
        type(COOc64_t), intent(in) :: COO
        complex(wp), intent(inout) :: dense(:,:)
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

        associate( nnz=>COO%nnz, num_rows=>COO%nrows, num_cols=>COO%ncols, base=>COO%base, sym=>COO%sym )
        CSR%NNZ = nnz; CSR%nrows = num_rows; CSR%ncols = num_cols
        CSR%base = base; CSR%sym = sym

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

    subroutine coo2csr_ordered_r_sp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COOr32_t), intent(in)    :: COO
        type(CSRr32_t), intent(inout) :: CSR
        integer :: i

        associate( nnz=>COO%nnz, num_rows=>COO%nrows, num_cols=>COO%ncols, base=>COO%base, sym=>COO%sym )
        CSR%NNZ = nnz; CSR%nrows = num_rows; CSR%ncols = num_cols
        CSR%base = base; CSR%sym = sym

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

    subroutine coo2csr_ordered_r_dp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COOr64_t), intent(in)    :: COO
        type(CSRr64_t), intent(inout) :: CSR
        integer :: i

        associate( nnz=>COO%nnz, num_rows=>COO%nrows, num_cols=>COO%ncols, base=>COO%base, sym=>COO%sym )
        CSR%NNZ = nnz; CSR%nrows = num_rows; CSR%ncols = num_cols
        CSR%base = base; CSR%sym = sym

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


    subroutine coo2csr_ordered_c_sp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COOc32_t), intent(in)    :: COO
        type(CSRc32_t), intent(inout) :: CSR
        integer :: i

        associate( nnz=>COO%nnz, num_rows=>COO%nrows, num_cols=>COO%ncols, base=>COO%base, sym=>COO%sym )
        CSR%NNZ = nnz; CSR%nrows = num_rows; CSR%ncols = num_cols
        CSR%base = base; CSR%sym = sym

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

    subroutine coo2csr_ordered_c_dp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COOc64_t), intent(in)    :: COO
        type(CSRc64_t), intent(inout) :: CSR
        integer :: i

        associate( nnz=>COO%nnz, num_rows=>COO%nrows, num_cols=>COO%ncols, base=>COO%base, sym=>COO%sym )
        CSR%NNZ = nnz; CSR%nrows = num_rows; CSR%ncols = num_cols
        CSR%base = base; CSR%sym = sym

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

    
    subroutine csr2coo(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSR_t), intent(in)    :: CSR
        type(COO_t), intent(inout) :: COO
        integer :: i, j

        COO%NNZ = CSR%NNZ; COO%nrows = CSR%nrows; COO%ncols = CSR%ncols
        COO%base = CSR%base; COO%sym = CSR%sym

        if( .not.allocated(COO%index) ) allocate( COO%index(2,CSR%NNZ) )
        
        do i = 1, CSR%nrows
            do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                COO%index(1:2,j) = [i,CSR%col(j)]
            end do
        end do
    end subroutine

    subroutine csr2coo_r_sp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSRr32_t), intent(in)    :: CSR
        type(COOr32_t), intent(inout) :: COO
        integer :: i, j

        COO%NNZ = CSR%NNZ; COO%nrows = CSR%nrows; COO%ncols = CSR%ncols
        COO%base = CSR%base; COO%sym = CSR%sym

        if( .not.allocated(COO%data) ) then
            allocate( COO%data(CSR%NNZ) , source = CSR%data(1:CSR%NNZ) )
        else 
            COO%data(1:CSR%NNZ) = CSR%data(1:CSR%NNZ)
        end if

        if( .not.allocated(COO%index) ) allocate( COO%index(2,CSR%NNZ) )
        
        do i = 1, CSR%nrows
            do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                COO%index(1:2,j) = [i,CSR%col(j)]
            end do
        end do
    end subroutine

    subroutine csr2coo_r_dp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSRr64_t), intent(in)    :: CSR
        type(COOr64_t), intent(inout) :: COO
        integer :: i, j

        COO%NNZ = CSR%NNZ; COO%nrows = CSR%nrows; COO%ncols = CSR%ncols
        COO%base = CSR%base; COO%sym = CSR%sym

        if( .not.allocated(COO%data) ) then
            allocate( COO%data(CSR%NNZ) , source = CSR%data(1:CSR%NNZ) )
        else 
            COO%data(1:CSR%NNZ) = CSR%data(1:CSR%NNZ)
        end if

        if( .not.allocated(COO%index) ) allocate( COO%index(2,CSR%NNZ) )

        do i = 1, CSR%nrows
            do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                COO%index(1:2,j) = [i,CSR%col(j)]
            end do
        end do
    end subroutine

    subroutine csr2coo_c_sp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSRc32_t), intent(in)    :: CSR
        type(COOc32_t), intent(inout) :: COO
        integer :: i, j

        COO%NNZ = CSR%NNZ; COO%nrows = CSR%nrows; COO%ncols = CSR%ncols
        COO%base = CSR%base; COO%sym = CSR%sym

        if( .not.allocated(COO%data) ) then
            allocate( COO%data(CSR%NNZ) , source = CSR%data(1:CSR%NNZ) )
        else 
            COO%data(1:CSR%NNZ) = CSR%data(1:CSR%NNZ)
        end if

        if( .not.allocated(COO%index) ) allocate( COO%index(2,CSR%NNZ) )
        
        do i = 1, CSR%nrows
            do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                COO%index(1:2,j) = [i,CSR%col(j)]
            end do
        end do
    end subroutine

    subroutine csr2coo_c_dp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSRc64_t), intent(in)    :: CSR
        type(COOc64_t), intent(inout) :: COO
        integer :: i, j

        COO%NNZ = CSR%NNZ; COO%nrows = CSR%nrows; COO%ncols = CSR%ncols
        COO%base = CSR%base; COO%sym = CSR%sym

        if( .not.allocated(COO%data) ) then
            allocate( COO%data(CSR%NNZ) , source = CSR%data(1:CSR%NNZ) )
        else 
            COO%data(1:CSR%NNZ) = CSR%data(1:CSR%NNZ)
        end if

        if( .not.allocated(COO%index) ) allocate( COO%index(2,CSR%NNZ) )

        do i = 1, CSR%nrows
            do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                COO%index(1:2,j) = [i,CSR%col(j)]
            end do
        end do
    end subroutine
    
end module fsparse_conversions
