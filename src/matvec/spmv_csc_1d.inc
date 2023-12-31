integer, intent(in)     :: nrows
integer, intent(in)     :: ncols
integer, intent(in)     :: NNZ
integer, intent(in)     :: row(NNZ)
integer, intent(in)     :: colptr(ncols+1)
real(wp), intent(in)    :: data(NNZ)
real(wp), intent(in)    :: vec_x(ncols)
real(wp), intent(inout) :: vec_y(nrows)
integer :: i, j, adr1, adr2

do concurrent(j=1:ncols)
    adr1 = colptr(j)
    adr2 = colptr(j+1)-1
    do i = adr1, adr2
        vec_y(row(i)) = vec_y(row(i)) + data(i) * vec_x(j)
    end do
end do