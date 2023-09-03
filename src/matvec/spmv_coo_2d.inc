integer, intent(in)     :: N, M, NNZ, dim
integer, intent(in)     :: index(2,NNZ)
real(wp), intent(in)    :: data(NNZ)
real(wp), intent(in)    :: vec_x(dim,M)
real(wp), intent(inout) :: vec_y(dim,N)
integer :: k, ik, jk

do concurrent (k = 1:NNZ)
    ik = index(1,k)
    jk = index(2,k)
    vec_y(1:dim,ik) = vec_y(1:dim,ik) + data(k) * vec_x(1:dim,jk)
end do