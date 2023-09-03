integer, intent(in)     :: N, M, NNZ, dim
integer, intent(in)     :: col(NNZ)
integer, intent(in)     :: rowptr(N+1)
real(wp), intent(in)    :: data(NNZ)
real(wp), intent(in)    :: vec_x(dim,M)
real(wp), intent(inout) :: vec_y(dim,N)
integer :: i, j, adr1, adr2

do concurrent(i=1:N)
    adr1 = rowptr(i)
    adr2 = rowptr(i+1)-1
    vec_y(1:dim,i) = 0._wp
    do j = adr1, adr2
        vec_y(1:dim,i) = vec_y(1:dim,i) + data(j) * vec_x(1:dim,col(j))
    end do
end do