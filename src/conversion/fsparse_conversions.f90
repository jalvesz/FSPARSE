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

    interface dense2diagonal
        module procedure dense2diagonal_sp
        module procedure dense2diagonal_dp
        module procedure dense2diagonal_csp
        module procedure dense2diagonal_cdp
    end interface

    interface coo2dense
        module procedure coo2dense_sp
        module procedure coo2dense_dp
        module procedure coo2dense_csp
        module procedure coo2dense_cdp
    end interface

    interface coo2diagonal
        module procedure coo2diagonal_sp
        module procedure coo2diagonal_dp
        module procedure coo2diagonal_csp
        module procedure coo2diagonal_cdp
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

    interface csr2diagonal
        module procedure csr2diagonal_sp
        module procedure csr2diagonal_dp
        module procedure csr2diagonal_csp
        module procedure csr2diagonal_cdp
    end interface

    interface csr2sellc
        module procedure csr2sellc_sp
        module procedure csr2sellc_dp
        module procedure csr2sellc_csp
        module procedure csr2sellc_cdp
    end interface

    interface csr_block_expansion
        module procedure csr_block_expansion_sp
        module procedure csr_block_expansion_dp
        module procedure csr_block_expansion_csp
        module procedure csr_block_expansion_cdp
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


    subroutine dense2diagonal_sp(dense,diagonal)
        real(sp), intent(in) :: dense(:,:)
        real(sp), intent(inout) :: diagonal(:)
        integer :: num_rows
        integer :: i

        num_rows = size(dense,dim=1)
        do concurrent(i = 1:num_rows)
            diagonal(i) = dense(i,i)
        end do
    end subroutine

    subroutine dense2diagonal_dp(dense,diagonal)
        real(dp), intent(in) :: dense(:,:)
        real(dp), intent(inout) :: diagonal(:)
        integer :: num_rows
        integer :: i

        num_rows = size(dense,dim=1)
        do concurrent(i = 1:num_rows)
            diagonal(i) = dense(i,i)
        end do
    end subroutine

    subroutine dense2diagonal_csp(dense,diagonal)
        complex(sp), intent(in) :: dense(:,:)
        complex(sp), intent(inout) :: diagonal(:)
        integer :: num_rows
        integer :: i

        num_rows = size(dense,dim=1)
        do concurrent(i = 1:num_rows)
            diagonal(i) = dense(i,i)
        end do
    end subroutine

    subroutine dense2diagonal_cdp(dense,diagonal)
        complex(dp), intent(in) :: dense(:,:)
        complex(dp), intent(inout) :: diagonal(:)
        integer :: num_rows
        integer :: i

        num_rows = size(dense,dim=1)
        do concurrent(i = 1:num_rows)
            diagonal(i) = dense(i,i)
        end do
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


    subroutine coo2diagonal_sp(COO,diagonal)
        integer, parameter :: wp = sp
        type(COO_sp), intent(in) :: COO
        real(sp), intent(inout) :: diagonal(:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            if(COO%index(1,idx)==COO%index(2,idx)) &
            & diagonal( COO%index(1,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2diagonal_dp(COO,diagonal)
        integer, parameter :: wp = dp
        type(COO_dp), intent(in) :: COO
        real(dp), intent(inout) :: diagonal(:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            if(COO%index(1,idx)==COO%index(2,idx)) &
            & diagonal( COO%index(1,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2diagonal_csp(COO,diagonal)
        integer, parameter :: wp = sp
        type(COO_csp), intent(in) :: COO
        complex(sp), intent(inout) :: diagonal(:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            if(COO%index(1,idx)==COO%index(2,idx)) &
            & diagonal( COO%index(1,idx) ) = COO%data(idx)
        end do
    end subroutine

    subroutine coo2diagonal_cdp(COO,diagonal)
        integer, parameter :: wp = dp
        type(COO_cdp), intent(in) :: COO
        complex(dp), intent(inout) :: diagonal(:)
        integer :: idx

        do concurrent(idx = 1:COO%NNZ)
            if(COO%index(1,idx)==COO%index(2,idx)) &
            & diagonal( COO%index(1,idx) ) = COO%data(idx)
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


    subroutine csr2sellc_sp(CSR,SELLC,chunk)
        !! csr2sellc: This function enables transfering data from a CSR matrix to a SELL-C matrix
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(CSR_sp), intent(in)      :: CSR
        type(SELLC_sp), intent(inout) :: SELLC
        integer, intent(in), optional :: chunk
        real(sp), parameter :: zero = zero_sp
        integer :: i, j, num_chunks

        if(present(chunk)) SELLC%chunk_size = chunk

        SELLC%nrows = CSR%nrows; SELLC%ncols = CSR%ncols
        SELLC%base  = CSR%base;  SELLC%sym = CSR%sym
        associate( nrows=>SELLC%nrows, ncols=>SELLC%ncols, nnz=>SELLC%nnz, &
        &         chunk_size=>SELLC%chunk_size     )
        !-------------------------------------------
        ! csr rowptr to SELL-C chunked rowptr
        num_chunks = (nrows + chunk_size - 1)/chunk_size
        allocate( SELLC%rowptr(num_chunks+1) )
        block
            integer :: cidx, rownnz, chunknnz
            SELLC%rowptr(1) = 1
            cidx = 1
            do i = 1, nrows, chunk_size
                chunknnz = 0
                ! Iterate over rows in a given chunk
                do j = i, min(i+chunk_size-1,nrows)
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j)
                    chunknnz = max(chunknnz,rownnz)
                end do
                SELLC%rowptr(cidx+1) = SELLC%rowptr(cidx) + chunknnz
                cidx = cidx + 1
            end do
            nnz = SELLC%rowptr(num_chunks+1) - 1
        end block
        !-------------------------------------------
        ! copy values and colum index
        allocate(SELLC%col(chunk_size,nnz), source = -1)
        allocate(SELLC%data(chunk_size,nnz), source = zero )
        block
            integer :: lb, ri, iaa, iab, rownnz
            do i = 1, num_chunks

                lb = SELLC%rowptr(i)

                ! Loop over rows of a chunk
                do j = (i-1)*chunk_size + 1, min(i*chunk_size,nrows)
    
                    ri = j - (i - 1)*chunk_size
                    
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j) - 1
                    iaa    = CSR%rowptr(j)
                    iab    = CSR%rowptr(j+1) - 1
                    
                    SELLC%col(ri,lb:lb+rownnz)  = CSR%col(iaa:iab)
                    SELLC%data(ri,lb:lb+rownnz) = CSR%data(iaa:iab)
                
                end do
            end do
         end block
        end associate
    end subroutine

    subroutine csr2sellc_dp(CSR,SELLC,chunk)
        !! csr2sellc: This function enables transfering data from a CSR matrix to a SELL-C matrix
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(CSR_dp), intent(in)      :: CSR
        type(SELLC_dp), intent(inout) :: SELLC
        integer, intent(in), optional :: chunk
        real(dp), parameter :: zero = zero_dp
        integer :: i, j, num_chunks

        if(present(chunk)) SELLC%chunk_size = chunk

        SELLC%nrows = CSR%nrows; SELLC%ncols = CSR%ncols
        SELLC%base  = CSR%base;  SELLC%sym = CSR%sym
        associate( nrows=>SELLC%nrows, ncols=>SELLC%ncols, nnz=>SELLC%nnz, &
        &         chunk_size=>SELLC%chunk_size     )
        !-------------------------------------------
        ! csr rowptr to SELL-C chunked rowptr
        num_chunks = (nrows + chunk_size - 1)/chunk_size
        allocate( SELLC%rowptr(num_chunks+1) )
        block
            integer :: cidx, rownnz, chunknnz
            SELLC%rowptr(1) = 1
            cidx = 1
            do i = 1, nrows, chunk_size
                chunknnz = 0
                ! Iterate over rows in a given chunk
                do j = i, min(i+chunk_size-1,nrows)
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j)
                    chunknnz = max(chunknnz,rownnz)
                end do
                SELLC%rowptr(cidx+1) = SELLC%rowptr(cidx) + chunknnz
                cidx = cidx + 1
            end do
            nnz = SELLC%rowptr(num_chunks+1) - 1
        end block
        !-------------------------------------------
        ! copy values and colum index
        allocate(SELLC%col(chunk_size,nnz), source = -1)
        allocate(SELLC%data(chunk_size,nnz), source = zero )
        block
            integer :: lb, ri, iaa, iab, rownnz
            do i = 1, num_chunks

                lb = SELLC%rowptr(i)

                ! Loop over rows of a chunk
                do j = (i-1)*chunk_size + 1, min(i*chunk_size,nrows)
    
                    ri = j - (i - 1)*chunk_size
                    
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j) - 1
                    iaa    = CSR%rowptr(j)
                    iab    = CSR%rowptr(j+1) - 1
                    
                    SELLC%col(ri,lb:lb+rownnz)  = CSR%col(iaa:iab)
                    SELLC%data(ri,lb:lb+rownnz) = CSR%data(iaa:iab)
                
                end do
            end do
         end block
        end associate
    end subroutine

    subroutine csr2sellc_csp(CSR,SELLC,chunk)
        !! csr2sellc: This function enables transfering data from a CSR matrix to a SELL-C matrix
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(CSR_csp), intent(in)      :: CSR
        type(SELLC_csp), intent(inout) :: SELLC
        integer, intent(in), optional :: chunk
        complex(sp), parameter :: zero = zero_csp
        integer :: i, j, num_chunks

        if(present(chunk)) SELLC%chunk_size = chunk

        SELLC%nrows = CSR%nrows; SELLC%ncols = CSR%ncols
        SELLC%base  = CSR%base;  SELLC%sym = CSR%sym
        associate( nrows=>SELLC%nrows, ncols=>SELLC%ncols, nnz=>SELLC%nnz, &
        &         chunk_size=>SELLC%chunk_size     )
        !-------------------------------------------
        ! csr rowptr to SELL-C chunked rowptr
        num_chunks = (nrows + chunk_size - 1)/chunk_size
        allocate( SELLC%rowptr(num_chunks+1) )
        block
            integer :: cidx, rownnz, chunknnz
            SELLC%rowptr(1) = 1
            cidx = 1
            do i = 1, nrows, chunk_size
                chunknnz = 0
                ! Iterate over rows in a given chunk
                do j = i, min(i+chunk_size-1,nrows)
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j)
                    chunknnz = max(chunknnz,rownnz)
                end do
                SELLC%rowptr(cidx+1) = SELLC%rowptr(cidx) + chunknnz
                cidx = cidx + 1
            end do
            nnz = SELLC%rowptr(num_chunks+1) - 1
        end block
        !-------------------------------------------
        ! copy values and colum index
        allocate(SELLC%col(chunk_size,nnz), source = -1)
        allocate(SELLC%data(chunk_size,nnz), source = zero )
        block
            integer :: lb, ri, iaa, iab, rownnz
            do i = 1, num_chunks

                lb = SELLC%rowptr(i)

                ! Loop over rows of a chunk
                do j = (i-1)*chunk_size + 1, min(i*chunk_size,nrows)
    
                    ri = j - (i - 1)*chunk_size
                    
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j) - 1
                    iaa    = CSR%rowptr(j)
                    iab    = CSR%rowptr(j+1) - 1
                    
                    SELLC%col(ri,lb:lb+rownnz)  = CSR%col(iaa:iab)
                    SELLC%data(ri,lb:lb+rownnz) = CSR%data(iaa:iab)
                
                end do
            end do
         end block
        end associate
    end subroutine

    subroutine csr2sellc_cdp(CSR,SELLC,chunk)
        !! csr2sellc: This function enables transfering data from a CSR matrix to a SELL-C matrix
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(CSR_cdp), intent(in)      :: CSR
        type(SELLC_cdp), intent(inout) :: SELLC
        integer, intent(in), optional :: chunk
        complex(dp), parameter :: zero = zero_cdp
        integer :: i, j, num_chunks

        if(present(chunk)) SELLC%chunk_size = chunk

        SELLC%nrows = CSR%nrows; SELLC%ncols = CSR%ncols
        SELLC%base  = CSR%base;  SELLC%sym = CSR%sym
        associate( nrows=>SELLC%nrows, ncols=>SELLC%ncols, nnz=>SELLC%nnz, &
        &         chunk_size=>SELLC%chunk_size     )
        !-------------------------------------------
        ! csr rowptr to SELL-C chunked rowptr
        num_chunks = (nrows + chunk_size - 1)/chunk_size
        allocate( SELLC%rowptr(num_chunks+1) )
        block
            integer :: cidx, rownnz, chunknnz
            SELLC%rowptr(1) = 1
            cidx = 1
            do i = 1, nrows, chunk_size
                chunknnz = 0
                ! Iterate over rows in a given chunk
                do j = i, min(i+chunk_size-1,nrows)
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j)
                    chunknnz = max(chunknnz,rownnz)
                end do
                SELLC%rowptr(cidx+1) = SELLC%rowptr(cidx) + chunknnz
                cidx = cidx + 1
            end do
            nnz = SELLC%rowptr(num_chunks+1) - 1
        end block
        !-------------------------------------------
        ! copy values and colum index
        allocate(SELLC%col(chunk_size,nnz), source = -1)
        allocate(SELLC%data(chunk_size,nnz), source = zero )
        block
            integer :: lb, ri, iaa, iab, rownnz
            do i = 1, num_chunks

                lb = SELLC%rowptr(i)

                ! Loop over rows of a chunk
                do j = (i-1)*chunk_size + 1, min(i*chunk_size,nrows)
    
                    ri = j - (i - 1)*chunk_size
                    
                    rownnz = CSR%rowptr(j+1) - CSR%rowptr(j) - 1
                    iaa    = CSR%rowptr(j)
                    iab    = CSR%rowptr(j+1) - 1
                    
                    SELLC%col(ri,lb:lb+rownnz)  = CSR%col(iaa:iab)
                    SELLC%data(ri,lb:lb+rownnz) = CSR%data(iaa:iab)
                
                end do
            end do
         end block
        end associate
    end subroutine


    subroutine csr2diagonal_sp(CSR,diagonal)
        integer, parameter :: wp = sp
        type(CSR_sp), intent(in) :: CSR
        real(sp), intent(inout) :: diagonal(:)
        integer :: i, j
        select case(CSR%sym)
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i+1)-1 )
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i) )
            end do
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                    if( CSR%col(j) == i ) then
                        diagonal(i) = CSR%data(j)
                        exit
                    end if
                end do
            end do
        end select
    end subroutine

    subroutine csr2diagonal_dp(CSR,diagonal)
        integer, parameter :: wp = dp
        type(CSR_dp), intent(in) :: CSR
        real(dp), intent(inout) :: diagonal(:)
        integer :: i, j
        select case(CSR%sym)
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i+1)-1 )
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i) )
            end do
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                    if( CSR%col(j) == i ) then
                        diagonal(i) = CSR%data(j)
                        exit
                    end if
                end do
            end do
        end select
    end subroutine

    subroutine csr2diagonal_csp(CSR,diagonal)
        integer, parameter :: wp = sp
        type(CSR_csp), intent(in) :: CSR
        complex(sp), intent(inout) :: diagonal(:)
        integer :: i, j
        select case(CSR%sym)
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i+1)-1 )
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i) )
            end do
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                    if( CSR%col(j) == i ) then
                        diagonal(i) = CSR%data(j)
                        exit
                    end if
                end do
            end do
        end select
    end subroutine

    subroutine csr2diagonal_cdp(CSR,diagonal)
        integer, parameter :: wp = dp
        type(CSR_cdp), intent(in) :: CSR
        complex(dp), intent(inout) :: diagonal(:)
        integer :: i, j
        select case(CSR%sym)
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i+1)-1 )
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                diagonal(i) = CSR%data( CSR%rowptr(i) )
            end do
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                    if( CSR%col(j) == i ) then
                        diagonal(i) = CSR%data(j)
                        exit
                    end if
                end do
            end do
        end select
    end subroutine


    subroutine csr_block_expansion_sp(CSR,num_dof)
        !! This function enables expanding an ordered CSR matrix to include contigous blocks
        !! Warning: this operation does not preserve the data buffer
        type(CSR_sp), intent(inout) :: CSR
        integer, intent(in) :: num_dof
        real(sp), parameter :: zero = zero_sp
        integer, allocatable :: rowptr_expn(:), col_expn(:)
        integer :: i, j, p, q, adr1, adr2, block_nnz

        select case(CSR%sym)
        case(k_NOSYMMETRY)
            block_nnz = num_dof ** 2
            CSR%NNZ = CSR%NNZ * block_nnz
        case(k_SYMTRISUP,k_SYMTRIINF)
            block_nnz = num_dof + num_dof * (num_dof-1) / 2
            CSR%NNZ = CSR%nrows * block_nnz + (CSR%NNZ-CSR%nrows) * num_dof ** 2
        end select

        allocate( col_expn(CSR%NNZ) )
        allocate( rowptr_expn(CSR%nrows*num_dof+1) )

        rowptr_expn(1) = 1
        select case(CSR%sym)
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i))
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i))
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - p + 1
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    j = CSR%rowptr(i)
                    adr2 = adr1 - p + 1
                    do q = p, num_dof
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                    do j = CSR%rowptr(i)+1, CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) - p + 1
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) &
                    & + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - num_dof + p
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-2
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                    j = CSR%rowptr(i+1)-1
                    adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                    do q = 1, p
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                end do
            end do
        end select

        CSR%nrows = CSR%nrows * num_dof
        CSR%ncols = CSR%ncols * num_dof

        call move_alloc( rowptr_expn , CSR%rowptr )
        call move_alloc( col_expn    , CSR%col    )

        if( allocated(CSR%data) ) deallocate( CSR%data )
        allocate( CSR%data(CSR%NNZ) , source = zero )
    end subroutine

    subroutine csr_block_expansion_dp(CSR,num_dof)
        !! This function enables expanding an ordered CSR matrix to include contigous blocks
        !! Warning: this operation does not preserve the data buffer
        type(CSR_dp), intent(inout) :: CSR
        integer, intent(in) :: num_dof
        real(dp), parameter :: zero = zero_dp
        integer, allocatable :: rowptr_expn(:), col_expn(:)
        integer :: i, j, p, q, adr1, adr2, block_nnz

        select case(CSR%sym)
        case(k_NOSYMMETRY)
            block_nnz = num_dof ** 2
            CSR%NNZ = CSR%NNZ * block_nnz
        case(k_SYMTRISUP,k_SYMTRIINF)
            block_nnz = num_dof + num_dof * (num_dof-1) / 2
            CSR%NNZ = CSR%nrows * block_nnz + (CSR%NNZ-CSR%nrows) * num_dof ** 2
        end select

        allocate( col_expn(CSR%NNZ) )
        allocate( rowptr_expn(CSR%nrows*num_dof+1) )

        rowptr_expn(1) = 1
        select case(CSR%sym)
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i))
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i))
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - p + 1
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    j = CSR%rowptr(i)
                    adr2 = adr1 - p + 1
                    do q = p, num_dof
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                    do j = CSR%rowptr(i)+1, CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) - p + 1
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) &
                    & + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - num_dof + p
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-2
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                    j = CSR%rowptr(i+1)-1
                    adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                    do q = 1, p
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                end do
            end do
        end select

        CSR%nrows = CSR%nrows * num_dof
        CSR%ncols = CSR%ncols * num_dof

        call move_alloc( rowptr_expn , CSR%rowptr )
        call move_alloc( col_expn    , CSR%col    )

        if( allocated(CSR%data) ) deallocate( CSR%data )
        allocate( CSR%data(CSR%NNZ) , source = zero )
    end subroutine

    subroutine csr_block_expansion_csp(CSR,num_dof)
        !! This function enables expanding an ordered CSR matrix to include contigous blocks
        !! Warning: this operation does not preserve the data buffer
        type(CSR_csp), intent(inout) :: CSR
        integer, intent(in) :: num_dof
        complex(sp), parameter :: zero = zero_csp
        integer, allocatable :: rowptr_expn(:), col_expn(:)
        integer :: i, j, p, q, adr1, adr2, block_nnz

        select case(CSR%sym)
        case(k_NOSYMMETRY)
            block_nnz = num_dof ** 2
            CSR%NNZ = CSR%NNZ * block_nnz
        case(k_SYMTRISUP,k_SYMTRIINF)
            block_nnz = num_dof + num_dof * (num_dof-1) / 2
            CSR%NNZ = CSR%nrows * block_nnz + (CSR%NNZ-CSR%nrows) * num_dof ** 2
        end select

        allocate( col_expn(CSR%NNZ) )
        allocate( rowptr_expn(CSR%nrows*num_dof+1) )

        rowptr_expn(1) = 1
        select case(CSR%sym)
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i))
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i))
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - p + 1
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    j = CSR%rowptr(i)
                    adr2 = adr1 - p + 1
                    do q = p, num_dof
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                    do j = CSR%rowptr(i)+1, CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) - p + 1
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) &
                    & + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - num_dof + p
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-2
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                    j = CSR%rowptr(i+1)-1
                    adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                    do q = 1, p
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                end do
            end do
        end select

        CSR%nrows = CSR%nrows * num_dof
        CSR%ncols = CSR%ncols * num_dof

        call move_alloc( rowptr_expn , CSR%rowptr )
        call move_alloc( col_expn    , CSR%col    )

        if( allocated(CSR%data) ) deallocate( CSR%data )
        allocate( CSR%data(CSR%NNZ) , source = zero )
    end subroutine

    subroutine csr_block_expansion_cdp(CSR,num_dof)
        !! This function enables expanding an ordered CSR matrix to include contigous blocks
        !! Warning: this operation does not preserve the data buffer
        type(CSR_cdp), intent(inout) :: CSR
        integer, intent(in) :: num_dof
        complex(dp), parameter :: zero = zero_cdp
        integer, allocatable :: rowptr_expn(:), col_expn(:)
        integer :: i, j, p, q, adr1, adr2, block_nnz

        select case(CSR%sym)
        case(k_NOSYMMETRY)
            block_nnz = num_dof ** 2
            CSR%NNZ = CSR%NNZ * block_nnz
        case(k_SYMTRISUP,k_SYMTRIINF)
            block_nnz = num_dof + num_dof * (num_dof-1) / 2
            CSR%NNZ = CSR%nrows * block_nnz + (CSR%NNZ-CSR%nrows) * num_dof ** 2
        end select

        allocate( col_expn(CSR%NNZ) )
        allocate( rowptr_expn(CSR%nrows*num_dof+1) )

        rowptr_expn(1) = 1
        select case(CSR%sym)
        case(k_NOSYMMETRY)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i))
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i))
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRISUP)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - p + 1
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    j = CSR%rowptr(i)
                    adr2 = adr1 - p + 1
                    do q = p, num_dof
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                    do j = CSR%rowptr(i)+1, CSR%rowptr(i+1)-1
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) - p + 1
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                end do
            end do
        case(k_SYMTRIINF)
            do i = 1, CSR%nrows
                do p = 1, num_dof
                    rowptr_expn(num_dof*(i-1)+p+1) = rowptr_expn(num_dof*(i-1)+p) &
                    & + num_dof*(CSR%rowptr(i+1)-CSR%rowptr(i)) - num_dof + p
                end do
                do p = 1, num_dof
                    adr1 = rowptr_expn(num_dof*(i-1)+p)
                    do j = CSR%rowptr(i), CSR%rowptr(i+1)-2
                        adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                        do q = 1, num_dof
                            col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                        end do
                    end do
                    j = CSR%rowptr(i+1)-1
                    adr2 = adr1 + num_dof*(j-CSR%rowptr(i)) 
                    do q = 1, p
                        col_expn(adr2+q-1) = num_dof*(CSR%col(j)-1)+q
                    end do
                end do
            end do
        end select

        CSR%nrows = CSR%nrows * num_dof
        CSR%ncols = CSR%ncols * num_dof

        call move_alloc( rowptr_expn , CSR%rowptr )
        call move_alloc( col_expn    , CSR%col    )

        if( allocated(CSR%data) ) deallocate( CSR%data )
        allocate( CSR%data(CSR%NNZ) , source = zero )
    end subroutine

    
end module fsparse_conversions