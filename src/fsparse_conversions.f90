!>---------------------------------------------------
!> Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!>
!> Use of this source code is governed by a MIT
!> license that can be found in the LICENSE.md file
!>---------------------------------------------------
module fsparse_conversions
    use fsparse_constants
    use fsparse_matrix_gallery
    implicit none
    
    interface dense2coo
        module procedure dense2coo_sp
        module procedure dense2coo_dp
        module procedure dense2coo_csp
        module procedure dense2coo_cdp
    end interface

    interface coo2dense
        module procedure coo2dense_sp
        module procedure coo2dense_dp
        module procedure coo2dense_csp
        module procedure coo2dense_cdp
    end interface

    interface coo2csr
        module procedure coo2csr_ordered
        module procedure coo2csr_ordered_sp
        module procedure coo2csr_ordered_dp
        module procedure coo2csr_ordered_csp
        module procedure coo2csr_ordered_cdp
    end interface

    interface csr2coo
        module procedure csr2coo
        module procedure csr2coo_sp
        module procedure csr2coo_dp
        module procedure csr2coo_csp
        module procedure csr2coo_cdp
    end interface

contains
    subroutine dense2coo_sp(dense,COO)
        integer, parameter :: wp = sp
        real(sp), intent(in) :: dense(:,:)
        type(COO_sp), intent(inout) :: COO
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
        integer, parameter :: wp = dp
        real(dp), intent(in) :: dense(:,:)
        type(COO_dp), intent(inout) :: COO
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

    subroutine dense2coo_csp(dense,COO)
        integer, parameter :: wp = sp
        complex(sp), intent(in) :: dense(:,:)
        type(COO_csp), intent(inout) :: COO
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

    subroutine dense2coo_cdp(dense,COO)
        integer, parameter :: wp = dp
        complex(dp), intent(in) :: dense(:,:)
        type(COO_cdp), intent(inout) :: COO
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
        integer, parameter :: wp = sp
        type(COO_sp), intent(in) :: COO
        real(sp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2dense_dp(COO,dense)
        integer, parameter :: wp = dp
        type(COO_dp), intent(in) :: COO
        real(dp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2dense_csp(COO,dense)
        integer, parameter :: wp = sp
        type(COO_csp), intent(in) :: COO
        complex(sp), intent(inout) :: dense(:,:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            dense( COO%index(1,idx) , COO%index(2,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2dense_cdp(COO,dense)
        integer, parameter :: wp = dp
        type(COO_cdp), intent(in) :: COO
        complex(dp), intent(inout) :: dense(:,:)
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

    subroutine coo2csr_ordered_sp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COO_sp), intent(in)    :: COO
        type(CSR_sp), intent(inout) :: CSR
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

    subroutine coo2csr_ordered_dp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COO_dp), intent(in)    :: COO
        type(CSR_dp), intent(inout) :: CSR
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

    subroutine coo2csr_ordered_csp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COO_csp), intent(in)    :: COO
        type(CSR_csp), intent(inout) :: CSR
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

    subroutine coo2csr_ordered_cdp(COO,CSR)
        !! coo2csr_ordered: This function enables transfering data from a COO matrix to a CSR matrix
        !! under the hypothesis that the COO is already ordered.
        type(COO_cdp), intent(in)    :: COO
        type(CSR_cdp), intent(inout) :: CSR
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

    subroutine csr2coo_sp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSR_sp), intent(in)    :: CSR
        type(COO_sp), intent(inout) :: COO
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

    subroutine csr2coo_dp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSR_dp), intent(in)    :: CSR
        type(COO_dp), intent(inout) :: COO
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

    subroutine csr2coo_csp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSR_csp), intent(in)    :: CSR
        type(COO_csp), intent(inout) :: COO
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

    subroutine csr2coo_cdp(CSR,COO)
        !! csr2coo: This function enables transfering data from a CSR matrix to a COO matrix
        type(CSR_cdp), intent(in)    :: CSR
        type(COO_cdp), intent(inout) :: COO
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