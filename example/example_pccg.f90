program example_pccg
    use fsparse
    implicit none
    integer :: i, unit
    type(COO_dp) :: COO
    type(CSR_dp) :: CSR
    real(8), allocatable :: X(:), B(:)
    logical, allocatable :: dirichlet(:)
    !=======================================
    ! read matrix data
    open(newunit=unit,file='./example/coo.dat',form='formatted')
    read(unit,*) COO%nrows, COO%ncols, COO%nnz
    call COO%malloc(COO%nrows, COO%ncols, COO%nnz)
    allocate( X(COO%nrows), B(COO%nrows), dirichlet(COO%nrows))
    do i =1, COO%nnz 
        read(unit,*) COO%index(1:2,i), COO%data(i)
    end do
    close(unit)

    open(newunit=unit,file='./example/load.dat',form='formatted')
    read(unit,*) 
    read(unit,*) B(1:COO%nrows)
    close(unit)

    open(newunit=unit,file='./example/dirichlet.dat',form='formatted')
    read(unit,*) 
    read(unit,*) dirichlet(1:COO%nrows)
    close(unit)
    !=======================================
    ! Set solver
    call coo2csr(COO, CSR)
    block
        use fsparse_krylov
        call cgsolve(CSR,b,x,dirichlet)
    end block
end program