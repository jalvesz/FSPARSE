integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: NNZ
integer, intent(in)     :: index(2,NNZ)
complex(wp), intent(in)    :: data(NNZ)
complex(wp), intent(in)    :: vec_x(ncols)
complex(wp), intent(inout) :: vec_y(nrows)
integer :: k, ik, jk

do concurrent (k = 1:NNZ)
    ik = index(1,k)
    jk = index(2,k)
    vec_y(ik) = vec_y(ik) + data(k) * vec_x(jk)
    if( ik==jk ) cycle
    vec_y(jk) = vec_y(jk) + data(k) * vec_x(ik)
end do