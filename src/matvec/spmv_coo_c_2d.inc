integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: NNZ, dim
integer, intent(in)     :: index(2,NNZ)
complex(wp), intent(in)    :: data(NNZ)
complex(wp), intent(in)    :: vec_x(dim,nrows)
complex(wp), intent(inout) :: vec_y(dim,ncols)
integer :: k, ik, jk

do concurrent (k = 1:NNZ)
    ik = index(1,k)
    jk = index(2,k)
    vec_y(1:dim,ik) = vec_y(1:dim,ik) + data(k) * vec_x(1:dim,jk)
end do