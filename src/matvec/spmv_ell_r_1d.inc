integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: MNZ_P_ROW
integer, intent(in)     :: index(nrows,MNZ_P_ROW)
real(wp), intent(in)    :: data(nrows,MNZ_P_ROW)
real(wp), intent(in)    :: vec_x(ncols)
real(wp), intent(inout) :: vec_y(nrows)
integer :: row, col, k

do concurrent (row = 1:nrows, k = 1:MNZ_P_ROW)
    col = index(row,k)
    if(col>0) vec_y(row) = vec_y(row) + data(row,k) * vec_x(col)
end do