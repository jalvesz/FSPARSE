integer, intent(in)     :: N, M, NNZ
integer, intent(in)     :: col(NNZ)
integer, intent(in)     :: rowptr(N+1)
real(wp), intent(in)    :: data(NNZ)
real(wp), intent(in)    :: vec_x(M)
real(wp), intent(inout) :: vec_y(N)
real(wp) :: aux, aux2
integer :: i, j, adr1, adr2

do i = 1 , N
    adr1 = rowptr(i)
    adr2 = rowptr(i+1)-1
    
    aux2 = vec_x(i)
    aux  = aux2 * data(adr1)
    do j = adr1+1, adr2
        aux = aux + data(j) * vec_x(col(j))
        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
    end do
    vec_y(i) = vec_y(i) + aux
end do