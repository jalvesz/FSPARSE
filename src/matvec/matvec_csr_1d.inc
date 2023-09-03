real(wp), intent(in)    :: vec_x(:)
real(wp), intent(inout) :: vec_y(:)

if( matrix%sym == k_NOSYMMETRY) then
    call spmv_csr(matrix%data, matrix%col, matrix%rowptr, matrix%N, matrix%M, matrix%NNZ, vec_x, vec_y)
else 
    print *, 'to be implemented'
endif