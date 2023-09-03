integer, intent(in)     :: N, M, NNZ
integer, intent(in)     :: col(NNZ)
integer, intent(in)     :: rowptr(N+1)
real(wp), intent(in)    :: data(NNZ)
real(wp), intent(in)    :: vec_x(M)
real(wp), intent(inout) :: vec_y(N)
integer :: i, j, adr1, adr2

do concurrent(i=1:N)
    adr1 = rowptr(i)
    adr2 = rowptr(i+1)-1
    vec_y(i) = 0._wp
    do j = adr1, adr2
        vec_y(i) = vec_y(i) + data(j) * vec_x(col(j))
    end do
end do