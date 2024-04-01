complex(wp), intent(in)    :: vec_x(:,:)
complex(wp), intent(inout) :: vec_y(:,:)
integer :: dim

dim = size(vec_x,dim=1)
if( matrix%sym == k_NOSYMMETRY) then
    call spmv_coo(matrix%data, matrix%index, matrix%nrows, matrix%ncols, matrix%NNZ, dim, vec_x, vec_y)
else 
    call spmv_coo_sym(matrix%data, matrix%index, matrix%nrows, matrix%ncols, matrix%NNZ, dim, vec_x, vec_y)
endif