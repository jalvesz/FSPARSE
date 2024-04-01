integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: NNZ
integer, intent(in)     :: col(NNZ)
integer, intent(in)     :: rowptr(nrows+1)
complex(wp), intent(in)    :: data(NNZ)
complex(wp), intent(in)    :: vec_x(ncols)
complex(wp), intent(inout) :: vec_y(nrows)
complex(wp) :: aux, aux2
integer :: i, j, adr1, adr2

do i = 1 , nrows
    adr1 = rowptr(i)
    adr2 = rowptr(i+1)-1
    
    aux  = 0._wp
    aux2 = vec_x(i)
    do j = adr1, adr2-1
        aux = aux + data(j) * vec_x(col(j))
        vec_y(col(j)) = vec_y(col(j)) + data(j) * aux2
    end do
    aux = aux + data(j) * aux2
    vec_y(i) = vec_y(i) + aux
end do