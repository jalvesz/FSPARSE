module test_fsparse
    use iso_fortran_env, only: sp=>real32, dp=>real64
    use fsparse
    implicit none
    
    contains
    
    subroutine test_coo()
        real(sp), allocatable :: dense(:,:)
        type(COOr32_t) :: COO
        real(sp), allocatable :: vec_x(:)
        real(sp), allocatable :: vec_y(:)

        allocate( dense(4,5) , source = &
                reshape([9.0,4.0,0.0,4.0, &
                         0.0,7.0,8.0,0.0, &
                         0.0,0.0,-1.0,5.0,&
                         0.0,0.0,8.0,6.0,&
                         -3.0,0.0,0.0,0.0],[4,5]) )
        
        call dense2coo( dense , COO )
        
        allocate( vec_x(5) , source = 1._sp )
        allocate( vec_y(4) , source = 0._sp )
        
        vec_y = matmul( dense, vec_x )
        print *, vec_y ! 6.00000000       11.0000000       15.0000000       15.0000000

        vec_y = 0._sp
        call matvec( COO , vec_x, vec_y )
        print *, vec_y ! 6.00000000       11.0000000       15.0000000       15.0000000
        
    end subroutine

    subroutine test_csr()
        type(CSRr64_t) :: CSR
        real(dp), allocatable :: vec_x(:)
        real(dp), allocatable :: vec_y(:)

        call CSR%malloc(4,5,10)
        CSR%data(:)   = dble([9,-3,4,7,8,-1,8,4,5,6])
        CSR%col(:)    = [1,5,1,2,2,3,4,1,3,4]
        CSR%rowptr(:) = [1,3,5,8,11]
        
        allocate( vec_x(5) , source = 1._dp )
        allocate( vec_y(4) , source = 0._dp )
        call matvec( CSR , vec_x , vec_y )
        print *, vec_y
        
    end subroutine

    subroutine test_csc()
        type(CSCr64_t) :: CSC
        real(dp), allocatable :: vec_x(:)
        real(dp), allocatable :: vec_y(:)

        call CSC%malloc(4,5,10)
        CSC%data(:)   = dble([9,4,4,7,8,-1,5,8,6,-3])
        CSC%row(:)    = [1,2,4,2,3,3,4,3,4,1]
        CSC%colptr(:) = [1,4,6,8,10,11]
        
        allocate( vec_x(5) , source = 1._dp )
        allocate( vec_y(4) , source = 0._dp )
        call matvec( CSC , vec_x , vec_y )
        print *, vec_y
        
    end subroutine

    subroutine test_ell()
        type(ELLr32_t) :: ELL
        real(sp), allocatable :: vec_x(:)
        real(sp), allocatable :: vec_y(:)

        call ELL%malloc(4,5,3)
        ELL%data(1,1:3)   = real([9,-3,0])
        ELL%data(2,1:3)   = real([4,7,0])
        ELL%data(3,1:3)   = real([8,-1,8])
        ELL%data(4,1:3)   = real([4,5,6])

        ELL%index(1,1:3) = [1,5,0]
        ELL%index(2,1:3) = [1,2,0]
        ELL%index(3,1:3) = [2,3,4]
        ELL%index(4,1:3) = [1,3,4]
        
        allocate( vec_x(5) , source = 1._sp )
        allocate( vec_y(4) , source = 0._sp )
        call matvec( ELL , vec_x , vec_y )
        print *, vec_y
        
    end subroutine

end module