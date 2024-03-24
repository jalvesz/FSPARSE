integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: NNZ, dim
integer, intent(in)     :: col(NNZ)
integer, intent(in)     :: rowptr(nrows+1)
complex(wp), intent(in)    :: data(NNZ)
complex(wp), intent(in)    :: vec_x(dim,nrows)
complex(wp), intent(inout) :: vec_y(dim,ncols)
integer :: i, j, adr1, adr2

do concurrent(i=1:nrows)
    adr1 = rowptr(i)
    adr2 = rowptr(i+1)-1
    vec_y(1:dim,i) = 0._wp
    do j = adr1, adr2
        vec_y(1:dim,i) = vec_y(1:dim,i) + data(j) * vec_x(1:dim,col(j))
    end do
end do