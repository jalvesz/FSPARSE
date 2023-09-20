integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: MNZ_P_ROW, dim
integer, intent(in)     :: index(nrows,MNZ_P_ROW)
real(wp), intent(in)    :: data(nrows,MNZ_P_ROW)
real(wp), intent(in)    :: vec_x(dim,ncols)
real(wp), intent(inout) :: vec_y(dim,nrows)
integer :: row, col, k

do concurrent (row = 1:nrows, k = 1:MNZ_P_ROW)
    col = index(row,k)
    if(col>0) vec_y(1:dim,row) = vec_y(1:dim,row) + data(row,k) * vec_x(1:dim,col)
end do