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

        module procedure matvec_coo_2d_sp
        module procedure matvec_csr_2d_sp
        module procedure matvec_csc_2d_sp
        module procedure matvec_ell_2d_sp

        module procedure matvec_coo_1d_dp
        module procedure matvec_csr_1d_dp
        module procedure matvec_csc_1d_dp
        module procedure matvec_ell_1d_dp

        module procedure matvec_coo_2d_dp
        module procedure matvec_csr_2d_dp
        module procedure matvec_csc_2d_dp
        module procedure matvec_ell_2d_dp

        module procedure matvec_coo_1d_qp
        module procedure matvec_csr_1d_qp
        module procedure matvec_csc_1d_qp
        module procedure matvec_ell_1d_qp

        module procedure matvec_coo_2d_qp
        module procedure matvec_csr_2d_qp
        module procedure matvec_csc_2d_qp
        module procedure matvec_ell_2d_qp

        module procedure matvec_coo_1d_csp
        module procedure matvec_csr_1d_csp
        module procedure matvec_csc_1d_csp
        module procedure matvec_ell_1d_csp

        module procedure matvec_coo_2d_csp
        module procedure matvec_csr_2d_csp
        module procedure matvec_csc_2d_csp
        module procedure matvec_ell_2d_csp

        module procedure matvec_coo_1d_cdp
        module procedure matvec_csr_1d_cdp
        module procedure matvec_csc_1d_cdp
        module procedure matvec_ell_1d_cdp

        module procedure matvec_coo_2d_cdp
        module procedure matvec_csr_2d_cdp
        module procedure matvec_csc_2d_cdp
        module procedure matvec_ell_2d_cdp

        module procedure matvec_coo_1d_cqp
        module procedure matvec_csr_1d_cqp
        module procedure matvec_csc_1d_cqp
        module procedure matvec_ell_1d_cqp

        module procedure matvec_coo_2d_cqp
        module procedure matvec_csr_2d_cqp
        module procedure matvec_csc_2d_cqp
        module procedure matvec_ell_2d_cqp

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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + data(k) * vec_x(:,ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_1d_qp(matrix,vec_x,vec_y)
        type(COO_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:)
        real(qp), intent(inout) :: vec_y(:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + data(k) * vec_x(ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_2d_qp(matrix,vec_x,vec_y)
        type(COO_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:,:)
        real(qp), intent(inout) :: vec_y(:,:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                    if( ik==jk ) cycle
                    vec_y(:,jk) = vec_y(:,jk) + data(k) * vec_x(:,ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_1d_cqp(matrix,vec_x,vec_y)
        type(COO_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:)
        complex(qp), intent(inout) :: vec_y(:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                end do

            else 
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
                    if( ik==jk ) cycle
                    vec_y(jk) = vec_y(jk) + data(k) * vec_x(ik)
                end do

            end if
        end associate
    end subroutine

    subroutine matvec_coo_2d_cqp(matrix,vec_x,vec_y)
        type(COO_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:,:)
        complex(qp), intent(inout) :: vec_y(:,:)
        integer :: k, ik, jk

        associate( data => matrix%data, index => matrix%index, sym => matrix%sym, nnz => matrix%nnz )
            if( sym == k_NOSYMMETRY) then
                do concurrent (k = 1:nnz)
                    ik = index(1,k)
                    jk = index(2,k)
                    vec_y(:,ik) = vec_y(:,ik) + data(k) * vec_x(:,jk)
                end do

            else 
                do concurrent (k = 1:nnz)
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
                do concurrent(i=1:nrows)
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
                do concurrent(i=1:nrows)
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
                do concurrent(i=1:nrows)
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
                do concurrent(i=1:nrows)
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
    
    subroutine matvec_csr_1d_qp(matrix,vec_x,vec_y)
        type(CSR_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:)
        real(qp), intent(inout) :: vec_y(:)
        integer :: i, j
        real(qp) :: aux

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(i=1:nrows)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_qp
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
    
    subroutine matvec_csr_2d_qp(matrix,vec_x,vec_y)
        type(CSR_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:,:)
        real(qp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        real(qp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(i=1:nrows)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,i) = vec_y(:,i) + data(j) * vec_x(:,col(j))
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_qp
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
                do concurrent(i=1:nrows)
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
                do concurrent(i=1:nrows)
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
                do concurrent(i=1:nrows)
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
                do concurrent(i=1:nrows)
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
    
    subroutine matvec_csr_1d_cqp(matrix,vec_x,vec_y)
        type(CSR_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:)
        complex(qp), intent(inout) :: vec_y(:)
        integer :: i, j
        complex(qp) :: aux

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(i=1:nrows)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_cqp
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
    
    subroutine matvec_csr_2d_cqp(matrix,vec_x,vec_y)
        type(CSR_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:,:)
        complex(qp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        complex(qp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, col => matrix%col, rowptr => matrix%rowptr, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(i=1:nrows)
                    do j = rowptr(i), rowptr(i+1)-1
                        vec_y(:,i) = vec_y(:,i) + data(j) * vec_x(:,col(j))
                    end do
                end do

            else if( sym == k_SYMTRIINF )then
                do i = 1 , nrows
                    aux  = zero_cqp
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
                do concurrent(j=1:ncols)
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
                do concurrent(j=1:ncols)
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
                do concurrent(j=1:ncols)
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
                do concurrent(j=1:ncols)
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
    
    subroutine matvec_csc_1d_qp(matrix,vec_x,vec_y)
        type(CSC_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:)
        real(qp), intent(inout) :: vec_y(:)
        integer :: i, j
        real(qp) :: aux

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(j=1:ncols)
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
                    aux  = zero_qp
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
    
    subroutine matvec_csc_2d_qp(matrix,vec_x,vec_y)
        type(CSC_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:,:)
        real(qp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        real(qp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(j=1:ncols)
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
                    aux  = zero_qp
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
                do concurrent(j=1:ncols)
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
                do concurrent(j=1:ncols)
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
                do concurrent(j=1:ncols)
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
                do concurrent(j=1:ncols)
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
    
    subroutine matvec_csc_1d_cqp(matrix,vec_x,vec_y)
        type(CSC_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:)
        complex(qp), intent(inout) :: vec_y(:)
        integer :: i, j
        complex(qp) :: aux

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(j=1:ncols)
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
                    aux  = zero_cqp
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
    
    subroutine matvec_csc_2d_cqp(matrix,vec_x,vec_y)
        type(CSC_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:,:)
        complex(qp), intent(inout) :: vec_y(:,:)
        integer :: i, j
        complex(qp) :: aux(size(vec_x,dim=1))

        associate( data => matrix%data, colptr => matrix%colptr, row => matrix%row, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent(j=1:ncols)
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
                    aux  = zero_cqp
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_1d_qp(matrix,vec_x,vec_y)
        type(ELL_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:)
        real(qp), intent(inout) :: vec_y(:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_2d_qp(matrix,vec_x,vec_y)
        type(ELL_qp), intent(in) :: matrix
        real(qp), intent(in)    :: vec_x(:,:)
        real(qp), intent(inout) :: vec_y(:,:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
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
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_1d_cqp(matrix,vec_x,vec_y)
        type(ELL_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:)
        complex(qp), intent(inout) :: vec_y(:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(i) = vec_y(i) + data(i,k) * vec_x(j)
                end do
                
            end if
        end associate
    end subroutine
    
    subroutine matvec_ell_2d_cqp(matrix,vec_x,vec_y)
        type(ELL_cqp), intent(in) :: matrix
        complex(qp), intent(in)    :: vec_x(:,:)
        complex(qp), intent(inout) :: vec_y(:,:)
        integer :: i, j, k

        associate( data => matrix%data, index => matrix%index, MNZ_P_ROW => matrix%K, &
            & nnz => matrix%nnz, nrows => matrix%nrows, ncols => matrix%ncols, sym => matrix%sym )
            if( sym == k_NOSYMMETRY) then
                do concurrent (i = 1:nrows, k = 1:MNZ_P_ROW)
                    j = index(i,k)
                    if(j>0) vec_y(:,i) = vec_y(:,i) + data(i,k) * vec_x(:,j)
                end do
                
            end if
        end associate
    end subroutine
    
    
end module fsparse_matvec