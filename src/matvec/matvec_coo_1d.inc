real(wp), intent(in)    :: vec_x(:)
real(wp), intent(inout) :: vec_y(:)

if( matrix%sym == k_NOSYMMETRY) then
    call spmv_coo(matrix%data, matrix%index, matrix%N, matrix%M, matrix%NNZ, vec_x, vec_y)
else 
    call spmv_coo_sym(matrix%data, matrix%index, matrix%N, matrix%M, matrix%NNZ, vec_x, vec_y)
endif