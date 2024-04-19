!---------------------------------------------------
! Copyright 2023-present Transvalor S.A.
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module fsparse_sort
    use fsparse_constants
    use fsparse_matrix_gallery
    implicit none
    private

    interface sort_coo
        module procedure sort_coo_unique
        module procedure sort_coo_unique_sp
        module procedure sort_coo_unique_dp
        module procedure sort_coo_unique_csp
        module procedure sort_coo_unique_cdp
    end interface

    interface quicksort
        module procedure quicksort_i
        module procedure quicksort_i_sp
        module procedure quicksort_i_dp
        module procedure quicksort_i_csp
        module procedure quicksort_i_cdp
    end interface

    public :: coo2ordered

    contains

    recursive subroutine quicksort_i(a, first, last)
        ! ref: https://gist.github.com/t-nissie/479f0f16966925fa29ea
        integer, intent(inout) :: a(*)
        integer, intent(in)    :: first, last
        integer :: i, j, x, t

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
                i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
            t = a(i);  a(i) = a(j);  a(j) = t
            i=i+1
            j=j-1
        end do
        if (first < i-1) call quicksort_i(a, first, i-1)
        if (j+1 < last)  call quicksort_i(a, j+1, last)
    end subroutine

    recursive subroutine quicksort_i_sp(a, b, first, last)
        integer, parameter :: wp = sp
        integer, intent(inout)  :: a(*) !! reference table to sort
        real(sp), intent(inout) :: b(*) !! secondary real data to sort w.r.t. a(:)
        integer, intent(in)     :: first, last
        integer  :: i, j, x, t
        real(sp) :: d

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
                i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
            t = a(i);  a(i) = a(j);  a(j) = t
            d = b(i);  b(i) = b(j);  b(j) = d
            i=i+1
            j=j-1
        end do
        if (first < i-1) call quicksort_i_sp(a, b, first, i-1)
        if (j+1 < last)  call quicksort_i_sp(a, b, j+1, last)
    end subroutine 

    recursive subroutine quicksort_i_dp(a, b, first, last)
        integer, parameter :: wp = sp
        integer, intent(inout)  :: a(*) !! reference table to sort
        real(dp), intent(inout) :: b(*) !! secondary real data to sort w.r.t. a(:)
        integer, intent(in)     :: first, last
        integer  :: i, j, x, t
        real(dp) :: d

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
                i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
            t = a(i);  a(i) = a(j);  a(j) = t
            d = b(i);  b(i) = b(j);  b(j) = d
            i=i+1
            j=j-1
        end do
        if (first < i-1) call quicksort_i_dp(a, b, first, i-1)
        if (j+1 < last)  call quicksort_i_dp(a, b, j+1, last)
    end subroutine 

    recursive subroutine quicksort_i_csp(a, b, first, last)
        integer, parameter :: wp = sp
        integer, intent(inout)  :: a(*) !! reference table to sort
        complex(sp), intent(inout) :: b(*) !! secondary real data to sort w.r.t. a(:)
        integer, intent(in)     :: first, last
        integer  :: i, j, x, t
        complex(sp) :: d

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
                i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
            t = a(i);  a(i) = a(j);  a(j) = t
            d = b(i);  b(i) = b(j);  b(j) = d
            i=i+1
            j=j-1
        end do
        if (first < i-1) call quicksort_i_csp(a, b, first, i-1)
        if (j+1 < last)  call quicksort_i_csp(a, b, j+1, last)
    end subroutine 

    recursive subroutine quicksort_i_cdp(a, b, first, last)
        integer, parameter :: wp = sp
        integer, intent(inout)  :: a(*) !! reference table to sort
        complex(dp), intent(inout) :: b(*) !! secondary real data to sort w.r.t. a(:)
        integer, intent(in)     :: first, last
        integer  :: i, j, x, t
        complex(dp) :: d

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
                i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
            t = a(i);  a(i) = a(j);  a(j) = t
            d = b(i);  b(i) = b(j);  b(j) = d
            i=i+1
            j=j-1
        end do
        if (first < i-1) call quicksort_i_cdp(a, b, first, i-1)
        if (j+1 < last)  call quicksort_i_cdp(a, b, j+1, last)
    end subroutine 


    subroutine sort_coo_unique( a, n, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), rows_(:), cols_(:)
        !---------------------------------------------------------
        ! Sort a first time with respect to first index using Count sort
        allocate( count_i( 0:num_rows ) , source = 0 )
        do ed = 1, n
            count_i( a(1,ed) ) = count_i( a(1,ed) ) + 1
        end do
        do n_i = 2, num_rows
            count_i(n_i) = count_i(n_i) + count_i(n_i-1)
        end do
        allocate( count_i_aux( 0:num_rows ) , source = count_i )

        allocate( rows_(n), cols_(n) )
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            rows_(pos) = a(1,ed)
            cols_(pos) = a(2,ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        !---------------------------------------------------------
        ! Sort with respect to second colum
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            if(dd>0) call quicksort_i(cols_(adr0),1,dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
        forall(ed=1:n) a(1:2,ed) = [rows_(ed),cols_(ed)]
        stride = 0
        do ed = 2, n
            if( a(1,ed) == a(1,ed-1) .and. a(2,ed) == a(2,ed-1) ) then
                stride = stride + 1
            else
                a(1:2,ed-stride) = a(1:2,ed)
            end if
        end do
        n = n - stride
    end subroutine

    subroutine sort_coo_unique_sp( a, data, n, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, parameter :: wp = sp
        real(sp), intent(inout) :: data(*)
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), rows_(:), cols_(:)
        real(sp), allocatable :: temp(:)
        !---------------------------------------------------------
        ! Sort a first time with respect to first index using Count sort
        allocate( count_i( 0:num_rows ) , source = 0 )
        do ed = 1, n
            count_i( a(1,ed) ) = count_i( a(1,ed) ) + 1
        end do
        do n_i = 2, num_rows
            count_i(n_i) = count_i(n_i) + count_i(n_i-1)
        end do
        allocate( count_i_aux( 0:num_rows ) , source = count_i )

        allocate( rows_(n), cols_(n), temp(n) )
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            rows_(pos) = a(1,ed)
            cols_(pos) = a(2,ed)
            temp(pos)  = data(ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        !---------------------------------------------------------
        ! Sort with respect to second colum using a quicksort
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            if(dd>0) call quicksort_i_sp(cols_(adr0),temp(adr0),1,dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
        do ed = 1,n
            a(1:2,ed) = [rows_(ed),cols_(ed)]
        end do
        data(1:n) = temp(1:n)
        stride = 0
        do ed = 2, n
            if( a(1,ed) == a(1,ed-1) .and. a(2,ed) == a(2,ed-1) ) then
                data(ed-1-stride) = data(ed-1-stride) + data(ed) ; data(ed) = data(ed-1-stride)
                stride = stride + 1
            else
                a(1:2,ed-stride) = a(1:2,ed)
                data(ed-stride) = data(ed)
            end if
        end do
        n = n - stride
    end subroutine

    subroutine sort_coo_unique_dp( a, data, n, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, parameter :: wp = dp
        real(dp), intent(inout) :: data(*)
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), rows_(:), cols_(:)
        real(dp), allocatable :: temp(:)
        !---------------------------------------------------------
        ! Sort a first time with respect to first index using Count sort
        allocate( count_i( 0:num_rows ) , source = 0 )
        do ed = 1, n
            count_i( a(1,ed) ) = count_i( a(1,ed) ) + 1
        end do
        do n_i = 2, num_rows
            count_i(n_i) = count_i(n_i) + count_i(n_i-1)
        end do
        allocate( count_i_aux( 0:num_rows ) , source = count_i )

        allocate( rows_(n), cols_(n), temp(n) )
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            rows_(pos) = a(1,ed)
            cols_(pos) = a(2,ed)
            temp(pos)  = data(ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        !---------------------------------------------------------
        ! Sort with respect to second colum using a quicksort
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            if(dd>0) call quicksort_i_dp(cols_(adr0),temp(adr0),1,dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
        do ed = 1,n
            a(1:2,ed) = [rows_(ed),cols_(ed)]
        end do
        data(1:n) = temp(1:n)
        stride = 0
        do ed = 2, n
            if( a(1,ed) == a(1,ed-1) .and. a(2,ed) == a(2,ed-1) ) then
                data(ed-1-stride) = data(ed-1-stride) + data(ed) ; data(ed) = data(ed-1-stride)
                stride = stride + 1
            else
                a(1:2,ed-stride) = a(1:2,ed)
                data(ed-stride) = data(ed)
            end if
        end do
        n = n - stride
    end subroutine

    subroutine sort_coo_unique_csp( a, data, n, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, parameter :: wp = sp
        complex(sp), intent(inout) :: data(*)
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), rows_(:), cols_(:)
        complex(sp), allocatable :: temp(:)
        !---------------------------------------------------------
        ! Sort a first time with respect to first index using Count sort
        allocate( count_i( 0:num_rows ) , source = 0 )
        do ed = 1, n
            count_i( a(1,ed) ) = count_i( a(1,ed) ) + 1
        end do
        do n_i = 2, num_rows
            count_i(n_i) = count_i(n_i) + count_i(n_i-1)
        end do
        allocate( count_i_aux( 0:num_rows ) , source = count_i )

        allocate( rows_(n), cols_(n), temp(n) )
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            rows_(pos) = a(1,ed)
            cols_(pos) = a(2,ed)
            temp(pos)  = data(ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        !---------------------------------------------------------
        ! Sort with respect to second colum using a quicksort
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            if(dd>0) call quicksort_i_csp(cols_(adr0),temp(adr0),1,dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
        do ed = 1,n
            a(1:2,ed) = [rows_(ed),cols_(ed)]
        end do
        data(1:n) = temp(1:n)
        stride = 0
        do ed = 2, n
            if( a(1,ed) == a(1,ed-1) .and. a(2,ed) == a(2,ed-1) ) then
                data(ed-1-stride) = data(ed-1-stride) + data(ed) ; data(ed) = data(ed-1-stride)
                stride = stride + 1
            else
                a(1:2,ed-stride) = a(1:2,ed)
                data(ed-stride) = data(ed)
            end if
        end do
        n = n - stride
    end subroutine

    subroutine sort_coo_unique_cdp( a, data, n, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, parameter :: wp = dp
        complex(dp), intent(inout) :: data(*)
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), rows_(:), cols_(:)
        complex(dp), allocatable :: temp(:)
        !---------------------------------------------------------
        ! Sort a first time with respect to first index using Count sort
        allocate( count_i( 0:num_rows ) , source = 0 )
        do ed = 1, n
            count_i( a(1,ed) ) = count_i( a(1,ed) ) + 1
        end do
        do n_i = 2, num_rows
            count_i(n_i) = count_i(n_i) + count_i(n_i-1)
        end do
        allocate( count_i_aux( 0:num_rows ) , source = count_i )

        allocate( rows_(n), cols_(n), temp(n) )
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            rows_(pos) = a(1,ed)
            cols_(pos) = a(2,ed)
            temp(pos)  = data(ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        !---------------------------------------------------------
        ! Sort with respect to second colum using a quicksort
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            if(dd>0) call quicksort_i_cdp(cols_(adr0),temp(adr0),1,dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
        do ed = 1,n
            a(1:2,ed) = [rows_(ed),cols_(ed)]
        end do
        data(1:n) = temp(1:n)
        stride = 0
        do ed = 2, n
            if( a(1,ed) == a(1,ed-1) .and. a(2,ed) == a(2,ed-1) ) then
                data(ed-1-stride) = data(ed-1-stride) + data(ed) ; data(ed) = data(ed-1-stride)
                stride = stride + 1
            else
                a(1:2,ed-stride) = a(1:2,ed)
                data(ed-stride) = data(ed)
            end if
        end do
        n = n - stride
    end subroutine


    subroutine coo2ordered(COO)
        class(COO_t), intent(inout) :: COO
        integer, allocatable :: itemp(:,:)
        
        if(COO%isOrdered) return
        
        select type (coo)
            type is( coo_t )
                call sort_coo(COO%index, COO%nnz, COO%nrows, COO%ncols)
            type is( coo_sp )
                block
                real(sp), allocatable :: temp(:)
                call sort_coo(COO%index, COO%data, COO%nnz, COO%nrows, COO%ncols)
                
                allocate( temp(COO%nnz) , source=COO%data(1:COO%nnz) )
                call move_alloc( temp , COO%data )
                end block
            type is( coo_dp )
                block
                real(dp), allocatable :: temp(:)
                call sort_coo(COO%index, COO%data, COO%nnz, COO%nrows, COO%ncols)
                
                allocate( temp(COO%nnz) , source=COO%data(1:COO%nnz) )
                call move_alloc( temp , COO%data )
                end block
            type is( coo_csp )
                block
                complex(sp), allocatable :: temp(:)
                call sort_coo(COO%index, COO%data, COO%nnz, COO%nrows, COO%ncols)
                
                allocate( temp(COO%nnz) , source=COO%data(1:COO%nnz) )
                call move_alloc( temp , COO%data )
                end block
            type is( coo_cdp )
                block
                complex(dp), allocatable :: temp(:)
                call sort_coo(COO%index, COO%data, COO%nnz, COO%nrows, COO%ncols)
                
                allocate( temp(COO%nnz) , source=COO%data(1:COO%nnz) )
                call move_alloc( temp , COO%data )
                end block
        end select
        
        allocate( itemp(2,COO%nnz) , source=COO%index(1:2,1:COO%nnz) )
        call move_alloc( itemp , COO%index )

        COO%isOrdered = .true.
    end subroutine
    
end module fsparse_sort