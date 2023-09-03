integer, intent(in)     :: N, M, MNZ_P_ROW, dim
integer, intent(in)     :: index(N,MNZ_P_ROW)
real(wp), intent(in)    :: data(N,MNZ_P_ROW)
real(wp), intent(in)    :: vec_x(dim,M)
real(wp), intent(inout) :: vec_y(dim,N)
integer :: row, col, k

do concurrent (row = 1:N, k = 1:MNZ_P_ROW)
    col = index(row,k)
    if(col>0) vec_y(1:dim,row) = vec_y(1:dim,row) + data(row,k) * vec_x(1:dim,col)
end do