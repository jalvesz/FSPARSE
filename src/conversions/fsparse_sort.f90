!---------------------------------------------------
! Copyright 2023-present Transvalor S.A.
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module fsparse_sort
    use iso_fortran_env
    use fsparse_matrix_gallery
    implicit none
    private

    interface sort_coo
        module procedure sort_coo_unique
        module procedure sort_coo_unique_sp
        module procedure sort_coo_unique_dp
    end interface

    interface coo2ordered
        module procedure coo2ordered_i
        module procedure coo2ordered_s
        module procedure coo2ordered_d
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
    end subroutine quicksort_i

    recursive subroutine quicksort_is(a, b, first, last)
        integer, parameter :: wp = real32
        integer, intent(inout)  :: a(*) !! reference table to sort
        real(wp), intent(inout) :: b(*) !! secondary real data to sort w.r.t. a(:)
        integer, intent(in)     :: first, last
        integer  :: i, j, x, t
        real(wp) :: d

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
        if (first < i-1) call quicksort_is(a, b, first, i-1)
        if (j+1 < last)  call quicksort_is(a, b, j+1, last)
    end subroutine quicksort_is

    recursive subroutine quicksort_id(a, b, first, last)
        integer, parameter :: wp = real64
        integer, intent(inout)  :: a(*) !! reference table to sort
        real(wp), intent(inout) :: b(*) !! secondary real data to sort w.r.t. a(:)
        integer, intent(in)     :: first, last
        integer  :: i, j, x, t
        real(wp) :: d

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
        if (first < i-1) call quicksort_id(a, b, first, i-1)
        if (j+1 < last)  call quicksort_id(a, b, j+1, last)
    end subroutine quicksort_id

    subroutine sort_coo_unique( a, n, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, maxi, maxj, sze, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), temp_a(:,:), map_j(:)
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

        allocate( temp_a( 2,n ) ) 
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            temp_a(1:2,pos) = a(1:2,ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        a(1:2,1:n) = temp_a(1:2,1:n)
        !---------------------------------------------------------
        ! Sort with respect to second colum using std::sort
        sze = num_cols
        if( num_cols < 100 ) sze = 2*num_cols
        allocate( map_j( sze ) , source = 0 )
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            map_j(1:dd) = a(2,adr0:adr1)
            call quicksort_i(map_j,1,dd)
            a(2,adr0:adr1) = map_j(1:dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
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

    subroutine sort_coo_unique_sp( a, n, data, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, parameter :: wp = real32 
        real(wp), intent(inout) :: data(*)
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, maxi, maxj, sze, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), temp_a(:,:), map_j(:)
        real(wp), allocatable :: rtemp(:)
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

        allocate( temp_a( 2,n ) , rtemp(n) ) 
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            temp_a(1:2,pos) = a(1:2,ed)
            rtemp(pos) = data(ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        a(1:2,1:n) = temp_a(1:2,1:n)
        data(1:n) = rtemp(1:n)
        !---------------------------------------------------------
        ! Sort with respect to second colum using a quicksort
        sze = num_cols
        if( num_cols < 100 ) sze = 2*num_cols
        allocate( map_j( sze ) , source = 0 )
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            map_j(1:dd) = a(2,adr0:adr1)
            call quicksort_is(map_j,data(adr0),1,dd)
            a(2,adr0:adr1) = map_j(1:dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
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

    subroutine sort_coo_unique_dp( a, n, data, num_rows, num_cols )
        !! Sort a 2d array in increasing order first by index 1 and then by index 2
        integer, parameter :: wp = real64
        real(wp), intent(inout) :: data(*)
        integer, intent(inout) :: a(2,*)
        integer, intent(inout) :: n
        integer, intent(in) :: num_rows
        integer, intent(in) :: num_cols

        integer :: stride, maxi, maxj, sze, adr0, adr1, dd
        integer :: n_i, pos, ed
        integer, allocatable :: count_i(:), count_i_aux(:), temp_a(:,:), map_j(:)
        real(wp), allocatable :: rtemp(:)
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

        allocate( temp_a( 2,n ) , rtemp(n) ) 
        do ed = n, 1, -1
            n_i = a(1,ed)
            pos = count_i(n_i)
            temp_a(1:2,pos) = a(1:2,ed)
            rtemp(pos) = data(ed)
            count_i(n_i) = count_i(n_i) - 1
        end do
        a(1:2,1:n) = temp_a(1:2,1:n)
        data(1:n) = rtemp(1:n)
        !---------------------------------------------------------
        ! Sort with respect to second colum using a quicksort
        sze = num_cols
        if( num_cols < 100 ) sze = 2*num_cols
        allocate( map_j( sze ) , source = 0 )
        do n_i = 1, num_rows
            adr0 = count_i_aux(n_i-1)+1
            adr1 = count_i_aux(n_i)
            dd = adr1-adr0+1
            map_j(1:dd) = a(2,adr0:adr1)
            call quicksort_id(map_j,data(adr0),1,dd)
            a(2,adr0:adr1) = map_j(1:dd)
        end do
        !---------------------------------------------------------
        ! Remove duplicates
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

    subroutine coo2ordered_i(COO)
        type(COO_t), intent(inout) :: COO
        integer, allocatable :: itemp(:,:)

        if(.not.COO%isOrdered) then
            call sort_coo(COO%index, COO%nnz, COO%nrows, COO%ncols)

            allocate( itemp(2,COO%nnz) , source=COO%index(1:2,1:COO%nnz) )
            call move_alloc( itemp , COO%index )

            COO%isOrdered = .true.
        end if
    end subroutine

    subroutine coo2ordered_s(COO)
        integer, parameter :: wp=real32
        type(COOr32_t), intent(inout) :: COO
        integer, allocatable :: itemp(:,:)
        real(wp), allocatable :: rtemp(:)

        if(.not.COO%isOrdered) then
            call sort_coo(COO%index, COO%nnz, COO%data, COO%nrows, COO%ncols)

            allocate( itemp(2,COO%nnz) , source=COO%index(1:2,1:COO%nnz) )
            call move_alloc( itemp , COO%index )

            allocate( rtemp(COO%nnz) , source=COO%data(1:COO%nnz) )
            call move_alloc( rtemp , COO%data )

            COO%isOrdered = .true.
        end if
    end subroutine

    subroutine coo2ordered_d(COO)
        integer, parameter :: wp=real64
        type(COOr64_t), intent(inout) :: COO
        integer, allocatable :: itemp(:,:)
        real(wp), allocatable :: rtemp(:)

        if(.not.COO%isOrdered) then
            call sort_coo(COO%index, COO%nnz, COO%data, COO%nrows, COO%ncols)

            allocate( itemp(2,COO%nnz) , source=COO%index(1:2,1:COO%nnz) )
            call move_alloc( itemp , COO%index )

            allocate( rtemp(COO%nnz) , source=COO%data(1:COO%nnz) )
            call move_alloc( rtemp , COO%data )

            COO%isOrdered = .true.
        end if
    end subroutine
end module fsparse_sort