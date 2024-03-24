real(wp), intent(in)    :: vec_x(:)
real(wp), intent(inout) :: vec_y(:)

call spmv_ell(matrix%data, matrix%index, matrix%nrows, matrix%ncols, matrix%K, vec_x, vec_y)