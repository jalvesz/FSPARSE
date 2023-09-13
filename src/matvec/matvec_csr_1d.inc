real(wp), intent(in)    :: vec_x(:)
real(wp), intent(inout) :: vec_y(:)

if( matrix%sym == k_NOSYMMETRY) then
    call spmv_csr(matrix%data, matrix%col, matrix%rowptr, matrix%nrows, matrix%ncols, matrix%NNZ, vec_x, vec_y)
else if( matrix%sym == k_SYMTRIINF ) then
    call spmv_csr_syminf(matrix%data, matrix%col, matrix%rowptr, matrix%nrows, matrix%ncols, matrix%NNZ, vec_x, vec_y)
else if( matrix%sym == k_SYMTRISUP ) then
    call spmv_csr_symsup(matrix%data, matrix%col, matrix%rowptr, matrix%nrows, matrix%ncols, matrix%NNZ, vec_x, vec_y)
endif