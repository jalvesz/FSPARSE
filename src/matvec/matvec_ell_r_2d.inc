real(wp), intent(in)    :: vec_x(:,:)
real(wp), intent(inout) :: vec_y(:,:)
integer :: dim

dim = size(vec_x,dim=1)
call spmv_ell(matrix%data, matrix%index, matrix%nrows, matrix%ncols, matrix%K,dim, vec_x, vec_y)