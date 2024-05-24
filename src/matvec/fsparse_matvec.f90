!>---------------------------------------------------
!> Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!>
!> Use of this source code is governed by a MIT
!> license that can be found in the LICENSE.md file
!>---------------------------------------------------
module fsparse_matvec
    use fsparse_constants
    use fsparse_matrix_gallery
    implicit none
    private

    public :: matvec
    interface matvec ! call matvec(matrix,vec_x,vec_y) => vec_y = matrix * vec_x
        module procedure matvec_coo_1d_sp
        module procedure matvec_csr_1d_sp
        module procedure matvec_csc_1d_sp
        module procedure matvec_ell_1d_sp
        module procedure matvec_dense_1d_sp
        module procedure matvec_diagonal_1d_sp
        module procedure matvec_coo_2d_sp
        module procedure matvec_csr_2d_sp
        module procedure matvec_csc_2d_sp
        module procedure matvec_ell_2d_sp
        module procedure matvec_dense_2d_sp
        module procedure matvec_diagonal_2d_sp
        module procedure matvec_sellc_sp
        module procedure matvec_coo_1d_dp
        module procedure matvec_csr_1d_dp
        module procedure matvec_csc_1d_dp
        module procedure matvec_ell_1d_dp
        module procedure matvec_dense_1d_dp
        module procedure matvec_diagonal_1d_dp
        module procedure matvec_coo_2d_dp
        module procedure matvec_csr_2d_dp
        module procedure matvec_csc_2d_dp
        module procedure matvec_ell_2d_dp
        module procedure matvec_dense_2d_dp
        module procedure matvec_diagonal_2d_dp
        module procedure matvec_sellc_dp
        module procedure matvec_coo_1d_csp
        module procedure matvec_csr_1d_csp
        module procedure matvec_csc_1d_csp
        module procedure matvec_ell_1d_csp
        module procedure matvec_dense_1d_csp
        module procedure matvec_diagonal_1d_csp
        module procedure matvec_coo_2d_csp
        module procedure matvec_csr_2d_csp
        module procedure matvec_csc_2d_csp
        module procedure matvec_ell_2d_csp
        module procedure matvec_dense_2d_csp
        module procedure matvec_diagonal_2d_csp
        module procedure matvec_sellc_csp
        module procedure matvec_coo_1d_cdp
        module procedure matvec_csr_1d_cdp
        module procedure matvec_csc_1d_cdp
        module procedure matvec_ell_1d_cdp
        module procedure matvec_dense_1d_cdp
        module procedure matvec_diagonal_1d_cdp
        module procedure matvec_coo_2d_cdp
        module procedure matvec_csr_2d_cdp
        module procedure matvec_csc_2d_cdp
        module procedure matvec_ell_2d_cdp
        module procedure matvec_dense_2d_cdp
        module procedure matvec_diagonal_2d_cdp
        module procedure matvec_sellc_cdp
    end interface

contains

    !! matvec_coo
    subroutine matvec_coo_1d_sp(matrix,vec_x,vec_y)
        type(COO_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + data(k) * vec_x(ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_2d_sp(matrix,vec_x,vec_y)
        type(COO_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + data(k) * vec_x(:,ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_1d_dp(matrix,vec_x,vec_y)
        type(COO_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + data(k) * vec_x(ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_2d_dp(matrix,vec_x,vec_y)
        type(COO_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + data(k) * vec_x(:,ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_1d_csp(matrix,vec_x,vec_y)
        type(COO_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + data(k) * vec_x(ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_2d_csp(matrix,vec_x,vec_y)
        type(COO_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + data(k) * vec_x(:,ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_1d_cdp(matrix,vec_x,vec_y)
        type(COO_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + data(k) * vec_x(ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_2d_cdp(matrix,vec_x,vec_y)
        type(COO_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do k = 1, nnz
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + data(k) * vec_x(:,ik)
                end do

            end if
        end associate
    end subroutine


    !! matvec_csr
    subroutine matvec_csr_1d_sp(matrix,vec_x,vec_y)
        type(CSR_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        integer :: i, j
        real(sp) :: aux
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_sp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    aux = aux + data(j) * vec_x(i)
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    vec_y(i) = vec_y(i) + aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine matvec_csr_2d_sp(matrix,vec_x,vec_y)
        type(CSR_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        real(sp) :: aux(size(vec_x,dim=1))
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,i) = vec_y(:,i) + data(j) * vec_x(:,col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_sp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    aux = aux + data(j) * vec_x(:,i)
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine matvec_csr_1d_dp(matrix,vec_x,vec_y)
        type(CSR_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        integer :: i, j
        real(dp) :: aux
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_dp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    aux = aux + data(j) * vec_x(i)
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    vec_y(i) = vec_y(i) + aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine matvec_csr_2d_dp(matrix,vec_x,vec_y)
        type(CSR_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        real(dp) :: aux(size(vec_x,dim=1))
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,i) = vec_y(:,i) + data(j) * vec_x(:,col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_dp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    aux = aux + data(j) * vec_x(:,i)
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine matvec_csr_1d_csp(matrix,vec_x,vec_y)
        type(CSR_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        integer :: i, j
        complex(sp) :: aux
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_csp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    aux = aux + data(j) * vec_x(i)
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    vec_y(i) = vec_y(i) + aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine matvec_csr_2d_csp(matrix,vec_x,vec_y)
        type(CSR_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        complex(sp) :: aux(size(vec_x,dim=1))
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,i) = vec_y(:,i) + data(j) * vec_x(:,col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_csp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    aux = aux + data(j) * vec_x(:,i)
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine matvec_csr_1d_cdp(matrix,vec_x,vec_y)
        type(CSR_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        integer :: i, j
        complex(dp) :: aux
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_cdp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    aux = aux + data(j) * vec_x(i)
                    vec_y(i) = vec_y(i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(col(j))
                        vec_y(col(j)) = vec_y(col(j)) + data(j) * vec_x(i)
                    end do
                    vec_y(i) = vec_y(i) + aux
                end do

            end if
        end associate
    end subroutine
    
    subroutine matvec_csr_2d_cdp(matrix,vec_x,vec_y)
        type(CSR_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        complex(dp) :: aux(size(vec_x,dim=1))
        
        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,i) = vec_y(:,i) + data(j) * vec_x(:,col(j))
                    end do
                end do
                
            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_cdp
                    do j = rowptr(i), rowptr(i+1)-2
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    aux = aux + data(j) * vec_x(:,i)
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                do i = 1 , nrows
                    aux  = vec_x(:,i) * data(rowptr(i))
                    do j = rowptr(i)+1, rowptr(i+1)-1
                        aux = aux + data(j) * vec_x(:,col(j))
                        vec_y(:,col(j)) = vec_y(:,col(j)) + data(j) * vec_x(:,i)
                    end do
                    vec_y(:,i) = vec_y(:,i) + aux
                end do

            end if
        end associate
    end subroutine
    

    !! matvec_csc
    subroutine matvec_csc_1d_sp(matrix,vec_x,vec_y)
        type(CSC_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        integer :: i, j
        real(sp) :: aux

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_sp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_csc_2d_sp(matrix,vec_x,vec_y)
        type(CSC_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        real(sp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_sp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_csc_1d_dp(matrix,vec_x,vec_y)
        type(CSC_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        integer :: i, j
        real(dp) :: aux

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_dp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_csc_2d_dp(matrix,vec_x,vec_y)
        type(CSC_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        real(dp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_dp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_csc_1d_csp(matrix,vec_x,vec_y)
        type(CSC_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        integer :: i, j
        complex(sp) :: aux

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_csp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_csc_2d_csp(matrix,vec_x,vec_y)
        type(CSC_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        complex(sp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_csp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_csc_1d_cdp(matrix,vec_x,vec_y)
        type(CSC_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        integer :: i, j
        complex(dp) :: aux

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(row(i))
                        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
                    end do
                    vec_y(j) = vec_y(j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_cdp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(j)
                        vec_y(i) = vec_y(i) + data(i) * vec_x(row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(j)
                    vec_y(j) = vec_y(j) + aux
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_csc_2d_cdp(matrix,vec_x,vec_y)
        type(CSC_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        complex(dp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do j = 1, ncols
                    do i = colptr(j), colptr(j+1)-1
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = vec_x(:,j) * data(colptr(j))
                    do i = colptr(j)+1, colptr(j+1)-1
                        aux = aux + data(i) * vec_x(:,row(i))
                        vec_y(:,row(i)) = vec_y(:,row(i)) + data(i) * vec_x(:,j)
                    end do
                    vec_y(:,j) = vec_y(:,j) + aux
                end do

            else if( sym == k_SYMTRISUP )then
                ! NOT TESTED
                do j = 1 , ncols
                    aux  = zero_cdp
                    do i = colptr(j), colptr(i+1)-2
                        aux = aux + data(i) * vec_x(:,j)
                        vec_y(:,i) = vec_y(:,i) + data(i) * vec_x(:,row(i))
                    end do
                    aux = aux + data(colptr(j)) * vec_x(:,j)
                    vec_y(:,j) = vec_y(:,j) + aux
                end do
                
            end if
        end associate
    end subroutine
    

    !! matvec_ell
    subroutine matvec_ell_1d_sp(matrix,vec_x,vec_y)
        type(ELL_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_2d_sp(matrix,vec_x,vec_y)
        type(ELL_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_1d_dp(matrix,vec_x,vec_y)
        type(ELL_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_2d_dp(matrix,vec_x,vec_y)
        type(ELL_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_1d_csp(matrix,vec_x,vec_y)
        type(ELL_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_2d_csp(matrix,vec_x,vec_y)
        type(ELL_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_1d_cdp(matrix,vec_x,vec_y)
        type(ELL_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_2d_cdp(matrix,vec_x,vec_y)
        type(ELL_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do i = 1, nrows
                    do k = 1, MNZ_P_ROW
                        j = index(i,k)
                        if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
                    end do
                end do
                
            end if
        end associate
    end subroutine
    

    !! matvec_sellc
    subroutine matvec_sellc_sp(matrix,vec_x,vec_y)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_sp), intent(in) :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        real(sp), parameter     :: zero = zero_sp
        integer :: i, nz, rowidx, num_chunks, rm

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym  )
        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( sym == k_NOSYMMETRY) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case default
                print *, "error: chunk size not supported."
                return
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_remainder(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if

        end if
        end associate

    contains
    pure subroutine chunk_kernel_4(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        real(sp), intent(in)  :: a(4,nz), x(*)
        integer, intent(in) :: ja(4,nz)
        real(sp), intent(out) :: y(4)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_8(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        real(sp), intent(in)  :: a(8,nz), x(*)
        integer, intent(in) :: ja(8,nz)
        real(sp), intent(out) :: y(8)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_16(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        real(sp), intent(in)  :: a(16,nz), x(*)
        integer, intent(in) :: ja(16,nz)
        real(sp), intent(out) :: y(16)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine

    pure subroutine chunk_kernel_remainder(nz,cs,rm,a,ja,x,y)
        integer, intent(in), value      :: nz, cs, rm
        real(sp), intent(in)  :: a(cs,nz), x(*)
        integer, intent(in) :: ja(cs,nz)
        real(sp), intent(out) :: y(rm)
        integer :: j
        do j = 1, nz
            where(ja(1:rm,j) > 0) y = y + a(1:rm,j) * x(ja(1:rm,j))               
        end do
    end subroutine  

    end subroutine
    
    subroutine matvec_sellc_dp(matrix,vec_x,vec_y)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_dp), intent(in) :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        real(dp), parameter     :: zero = zero_dp
        integer :: i, nz, rowidx, num_chunks, rm

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym  )
        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( sym == k_NOSYMMETRY) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case default
                print *, "error: chunk size not supported."
                return
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_remainder(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if

        end if
        end associate

    contains
    pure subroutine chunk_kernel_4(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        real(dp), intent(in)  :: a(4,nz), x(*)
        integer, intent(in) :: ja(4,nz)
        real(dp), intent(out) :: y(4)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_8(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        real(dp), intent(in)  :: a(8,nz), x(*)
        integer, intent(in) :: ja(8,nz)
        real(dp), intent(out) :: y(8)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_16(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        real(dp), intent(in)  :: a(16,nz), x(*)
        integer, intent(in) :: ja(16,nz)
        real(dp), intent(out) :: y(16)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine

    pure subroutine chunk_kernel_remainder(nz,cs,rm,a,ja,x,y)
        integer, intent(in), value      :: nz, cs, rm
        real(dp), intent(in)  :: a(cs,nz), x(*)
        integer, intent(in) :: ja(cs,nz)
        real(dp), intent(out) :: y(rm)
        integer :: j
        do j = 1, nz
            where(ja(1:rm,j) > 0) y = y + a(1:rm,j) * x(ja(1:rm,j))               
        end do
    end subroutine  

    end subroutine
    
    subroutine matvec_sellc_csp(matrix,vec_x,vec_y)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_csp), intent(in) :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        complex(sp), parameter     :: zero = zero_csp
        integer :: i, nz, rowidx, num_chunks, rm

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym  )
        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( sym == k_NOSYMMETRY) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case default
                print *, "error: chunk size not supported."
                return
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_remainder(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if

        end if
        end associate

    contains
    pure subroutine chunk_kernel_4(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        complex(sp), intent(in)  :: a(4,nz), x(*)
        integer, intent(in) :: ja(4,nz)
        complex(sp), intent(out) :: y(4)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_8(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        complex(sp), intent(in)  :: a(8,nz), x(*)
        integer, intent(in) :: ja(8,nz)
        complex(sp), intent(out) :: y(8)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_16(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        complex(sp), intent(in)  :: a(16,nz), x(*)
        integer, intent(in) :: ja(16,nz)
        complex(sp), intent(out) :: y(16)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine

    pure subroutine chunk_kernel_remainder(nz,cs,rm,a,ja,x,y)
        integer, intent(in), value      :: nz, cs, rm
        complex(sp), intent(in)  :: a(cs,nz), x(*)
        integer, intent(in) :: ja(cs,nz)
        complex(sp), intent(out) :: y(rm)
        integer :: j
        do j = 1, nz
            where(ja(1:rm,j) > 0) y = y + a(1:rm,j) * x(ja(1:rm,j))               
        end do
    end subroutine  

    end subroutine
    
    subroutine matvec_sellc_cdp(matrix,vec_x,vec_y)
        !! This algorithm was gracefully provided by Ivan Privec and adapted by Jose Alves
        type(SELLC_cdp), intent(in) :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        complex(dp), parameter     :: zero = zero_cdp
        integer :: i, nz, rowidx, num_chunks, rm

        associate( data => matrix%data, ia => matrix%rowptr , ja => matrix%col, cs => matrix%chunk_size, &
        &   nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym  )
        num_chunks = nrows / cs
        rm = nrows - num_chunks * cs
        if( sym == k_NOSYMMETRY) then

            select case(cs)
            case(4)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*4 + 1 
                    call chunk_kernel_4(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(8)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*8 + 1 
                    call chunk_kernel_8(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case(16)
                do i = 1, num_chunks
                    nz = ia(i+1) - ia(i)
                    rowidx = (i - 1)*16 + 1 
                    call chunk_kernel_16(nz,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
                end do
            case default
                print *, "error: chunk size not supported."
                return
            end select
            
            ! remainder
            if(rm>0)then 
                i = num_chunks + 1 
                nz = ia(i+1) - ia(i)
                rowidx = (i - 1)*cs + 1
                call chunk_kernel_remainder(nz,cs,rm,data(:,ia(i)),ja(:,ia(i)),vec_x,vec_y(rowidx:))
            end if

        end if
        end associate

    contains
    pure subroutine chunk_kernel_4(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        complex(dp), intent(in)  :: a(4,nz), x(*)
        integer, intent(in) :: ja(4,nz)
        complex(dp), intent(out) :: y(4)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_8(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        complex(dp), intent(in)  :: a(8,nz), x(*)
        integer, intent(in) :: ja(8,nz)
        complex(dp), intent(out) :: y(8)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine
    pure subroutine chunk_kernel_16(nz,a,ja,x,y)
        integer, intent(in), value      :: nz
        complex(dp), intent(in)  :: a(16,nz), x(*)
        integer, intent(in) :: ja(16,nz)
        complex(dp), intent(out) :: y(16)
        integer :: j
        do j = 1, nz
            where(ja(:,j) > 0) y = y + a(:,j) * x(ja(:,j))               
        end do
    end subroutine

    pure subroutine chunk_kernel_remainder(nz,cs,rm,a,ja,x,y)
        integer, intent(in), value      :: nz, cs, rm
        complex(dp), intent(in)  :: a(cs,nz), x(*)
        integer, intent(in) :: ja(cs,nz)
        complex(dp), intent(out) :: y(rm)
        integer :: j
        do j = 1, nz
            where(ja(1:rm,j) > 0) y = y + a(1:rm,j) * x(ja(1:rm,j))               
        end do
    end subroutine  

    end subroutine
    

    !! matvec_dense
    subroutine matvec_dense_1d_sp(matrix,vec_x,vec_y)
        type(dense_sp), intent(in)    :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matmul(matrix%data,vec_x)
    end subroutine
    
    subroutine matvec_dense_2d_sp(matrix,vec_x,vec_y)
        type(dense_sp), intent(in)    :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=1)
            vec_y(i,:) = matmul(matrix%data,vec_x(i,:))
        end do
    end subroutine
    
    subroutine matvec_dense_1d_dp(matrix,vec_x,vec_y)
        type(dense_dp), intent(in)    :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matmul(matrix%data,vec_x)
    end subroutine
    
    subroutine matvec_dense_2d_dp(matrix,vec_x,vec_y)
        type(dense_dp), intent(in)    :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=1)
            vec_y(i,:) = matmul(matrix%data,vec_x(i,:))
        end do
    end subroutine
    
    subroutine matvec_dense_1d_csp(matrix,vec_x,vec_y)
        type(dense_csp), intent(in)    :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matmul(matrix%data,vec_x)
    end subroutine
    
    subroutine matvec_dense_2d_csp(matrix,vec_x,vec_y)
        type(dense_csp), intent(in)    :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=1)
            vec_y(i,:) = matmul(matrix%data,vec_x(i,:))
        end do
    end subroutine
    
    subroutine matvec_dense_1d_cdp(matrix,vec_x,vec_y)
        type(dense_cdp), intent(in)    :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matmul(matrix%data,vec_x)
    end subroutine
    
    subroutine matvec_dense_2d_cdp(matrix,vec_x,vec_y)
        type(dense_cdp), intent(in)    :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=1)
            vec_y(i,:) = matmul(matrix%data,vec_x(i,:))
        end do
    end subroutine
    

    !! matvec_diagonal
    subroutine matvec_diagonal_1d_sp(matrix,vec_x,vec_y)
        type(diagonal_sp), intent(in)    :: matrix
        real(sp), intent(in)    :: vec_x(:)
        real(sp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matrix%data * vec_x
    end subroutine
    
    subroutine matvec_diagonal_2d_sp(matrix,vec_x,vec_y)
        type(diagonal_sp), intent(in)    :: matrix
        real(sp), intent(in)    :: vec_x(:,:)
        real(sp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=2)
            vec_y(:,i) = matrix%data(i) * vec_x(:,i)
        end do
    end subroutine
    
    subroutine matvec_diagonal_1d_dp(matrix,vec_x,vec_y)
        type(diagonal_dp), intent(in)    :: matrix
        real(dp), intent(in)    :: vec_x(:)
        real(dp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matrix%data * vec_x
    end subroutine
    
    subroutine matvec_diagonal_2d_dp(matrix,vec_x,vec_y)
        type(diagonal_dp), intent(in)    :: matrix
        real(dp), intent(in)    :: vec_x(:,:)
        real(dp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=2)
            vec_y(:,i) = matrix%data(i) * vec_x(:,i)
        end do
    end subroutine
    
    subroutine matvec_diagonal_1d_csp(matrix,vec_x,vec_y)
        type(diagonal_csp), intent(in)    :: matrix
        complex(sp), intent(in)    :: vec_x(:)
        complex(sp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matrix%data * vec_x
    end subroutine
    
    subroutine matvec_diagonal_2d_csp(matrix,vec_x,vec_y)
        type(diagonal_csp), intent(in)    :: matrix
        complex(sp), intent(in)    :: vec_x(:,:)
        complex(sp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=2)
            vec_y(:,i) = matrix%data(i) * vec_x(:,i)
        end do
    end subroutine
    
    subroutine matvec_diagonal_1d_cdp(matrix,vec_x,vec_y)
        type(diagonal_cdp), intent(in)    :: matrix
        complex(dp), intent(in)    :: vec_x(:)
        complex(dp), intent(inout) :: vec_y(:)
        integer :: i
        vec_y = matrix%data * vec_x
    end subroutine
    
    subroutine matvec_diagonal_2d_cdp(matrix,vec_x,vec_y)
        type(diagonal_cdp), intent(in)    :: matrix
        complex(dp), intent(in)    :: vec_x(:,:)
        complex(dp), intent(inout) :: vec_y(:,:)
        integer :: i
        do i = 1, size(vec_x,dim=2)
            vec_y(:,i) = matrix%data(i) * vec_x(:,i)
        end do
    end subroutine
    
    
end module fsparse_matvec