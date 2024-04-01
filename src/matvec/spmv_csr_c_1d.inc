integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: NNZ
integer, intent(in)     :: col(NNZ)
integer, intent(in)     :: rowptr(nrows+1)
complex(wp), intent(in)    :: data(NNZ)
complex(wp), intent(in)    :: vec_x(ncols)
complex(wp), intent(inout) :: vec_y(nrows)
integer :: i, j, adr1, adr2

do concurrent(i=1:nrows)
    adr1 = rowptr(i)
    adr2 = rowptr(i+1)-1
    vec_y(i) = 0._wp
    do j = adr1, adr2
        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
    end do
end do