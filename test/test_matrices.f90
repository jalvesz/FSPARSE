!>---------------------------------------------------
!> Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!>
!> Use of this source code is governed by a MIT
!> license that can be found in the LICENSE.md file
!>---------------------------------------------------
module test_fsparse
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fsparse
    implicit none
    
    contains

    !> Collect all exported unit tests
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest('coo', test_coo), &
            new_unittest('coo2ordered', test_coo2ordered), &
            new_unittest('csr', test_csr), &
            new_unittest('csc', test_csc), &
            new_unittest('ell', test_ell),  &
            new_unittest('symmetries', test_symmetries), &
            new_unittest('cells2sparse', test_cells2sparse) &
        ]
    end subroutine

    subroutine test_coo(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            type(COO_sp) :: COO
            real(sp), allocatable :: dense(:,:)
            real(sp), allocatable :: vec_x(:)
            real(sp), allocatable :: vec_y1(:), vec_y2(:)

            allocate( dense(4,5) , source = &
                    reshape(real([9,4, 0,4, &
                                  0,7, 8,0, &
                                  0,0,-1,5, &
                                  0,0, 8,6, &
                                 -3,0, 0,0],kind=wp),[4,5]) )
            
            call dense2coo( dense , COO )
            
            allocate( vec_x(5)  , source = 1._wp )
            allocate( vec_y1(4) , source = 0._wp )
            allocate( vec_y2(4) , source = 0._wp )

            vec_y1 = matmul( dense, vec_x )
            
            call check(error, all(vec_y1 == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return

            call matvec( COO, vec_x, vec_y2 )
            call check(error, all(vec_y1 == vec_y2) )
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            type(COO_dp) :: COO
            real(dp), allocatable :: dense(:,:)
            real(dp), allocatable :: vec_x(:)
            real(dp), allocatable :: vec_y1(:), vec_y2(:)

            allocate( dense(4,5) , source = &
                    reshape(real([9,4, 0,4, &
                                  0,7, 8,0, &
                                  0,0,-1,5, &
                                  0,0, 8,6, &
                                 -3,0, 0,0],kind=wp),[4,5]) )
            
            call dense2coo( dense , COO )
            
            allocate( vec_x(5)  , source = 1._wp )
            allocate( vec_y1(4) , source = 0._wp )
            allocate( vec_y2(4) , source = 0._wp )

            vec_y1 = matmul( dense, vec_x )
            
            call check(error, all(vec_y1 == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return

            call matvec( COO, vec_x, vec_y2 )
            call check(error, all(vec_y1 == vec_y2) )
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_coo2ordered(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(COO_sp) :: COO

        call COO%malloc(4,4,12)
        COO%data(:) = 1
        COO%index(:,1) = [1,2]
        COO%index(:,2) = [1,3]
        COO%index(:,3) = [1,4]
        COO%index(:,4) = [2,3]
        COO%index(:,5) = [2,4]
        COO%index(:,6) = [3,4]

        COO%index(:,7) = [2,3]
        COO%index(:,8) = [2,4]
        COO%index(:,9) = [2,5]
        COO%index(:,10) = [3,4]
        COO%index(:,11) = [3,5]
        COO%index(:,12) = [4,5]

        call coo2ordered(COO)
        call check(error, COO%NNZ < 12 .and. COO%NNZ == 9 )
        if (allocated(error)) return

        call check(error, all(COO%data==[1,1,1,2,2,1,2,1,1])  )
        if (allocated(error)) return

        call check(error, all(COO%index(1,:)==[1,1,1,2,2,2,3,3,4])  )
        if (allocated(error)) return

        call check(error, all(COO%index(2,:)==[2,3,4,3,4,5,4,5,5])  )
        if (allocated(error)) return

    end subroutine 

    subroutine test_csr(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            type(CSR_sp) :: CSR
            real(sp), allocatable :: vec_x(:)
            real(sp), allocatable :: vec_y(:)

            call CSR%malloc(4,5,10)
            CSR%data(:)   = real([9,-3,4,7,8,-1,8,4,5,6],kind=wp)
            CSR%col(:)    = [1,5,1,2,2,3,4,1,3,4]
            CSR%rowptr(:) = [1,3,5,8,11]
            
            allocate( vec_x(5) , source = 1._wp )
            allocate( vec_y(4) , source = 0._wp )
            call matvec( CSR, vec_x, vec_y )
            
            call check(error, all(vec_y == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            type(CSR_dp) :: CSR
            real(dp), allocatable :: vec_x(:)
            real(dp), allocatable :: vec_y(:)

            call CSR%malloc(4,5,10)
            CSR%data(:)   = real([9,-3,4,7,8,-1,8,4,5,6],kind=wp)
            CSR%col(:)    = [1,5,1,2,2,3,4,1,3,4]
            CSR%rowptr(:) = [1,3,5,8,11]
            
            allocate( vec_x(5) , source = 1._wp )
            allocate( vec_y(4) , source = 0._wp )
            call matvec( CSR, vec_x, vec_y )
            
            call check(error, all(vec_y == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_csc(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            type(CSC_sp) :: CSC
            real(sp), allocatable :: vec_x(:)
            real(sp), allocatable :: vec_y(:)

            call CSC%malloc(4,5,10)
            CSC%data(:)   = real([9,4,4,7,8,-1,5,8,6,-3],kind=wp)
            CSC%row(:)    = [1,2,4,2,3,3,4,3,4,1]
            CSC%colptr(:) = [1,4,6,8,10,11]
            
            allocate( vec_x(5) , source = 1._wp )
            allocate( vec_y(4) , source = 0._wp )
            call matvec( CSC, vec_x, vec_y )
            
            call check(error, all(vec_y == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            type(CSC_dp) :: CSC
            real(dp), allocatable :: vec_x(:)
            real(dp), allocatable :: vec_y(:)

            call CSC%malloc(4,5,10)
            CSC%data(:)   = real([9,4,4,7,8,-1,5,8,6,-3],kind=wp)
            CSC%row(:)    = [1,2,4,2,3,3,4,3,4,1]
            CSC%colptr(:) = [1,4,6,8,10,11]
            
            allocate( vec_x(5) , source = 1._wp )
            allocate( vec_y(4) , source = 0._wp )
            call matvec( CSC, vec_x, vec_y )
            
            call check(error, all(vec_y == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_ell(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            type(ELL_sp) :: ELL
            real(sp), allocatable :: vec_x(:)
            real(sp), allocatable :: vec_y(:)

            call ELL%malloc(4,5,10)
            ELL%data(1,1:3)   = real([9,-3,0],kind=wp)
            ELL%data(2,1:3)   = real([4,7,0],kind=wp)
            ELL%data(3,1:3)   = real([8,-1,8],kind=wp)
            ELL%data(4,1:3)   = real([4,5,6],kind=wp)

            ELL%index(1,1:3) = [1,5,0]
            ELL%index(2,1:3) = [1,2,0]
            ELL%index(3,1:3) = [2,3,4]
            ELL%index(4,1:3) = [1,3,4]
            
            allocate( vec_x(5) , source = 1._wp )
            allocate( vec_y(4) , source = 0._wp )
            call matvec( ELL, vec_x, vec_y )
            
            call check(error, all(vec_y == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            type(ELL_dp) :: ELL
            real(dp), allocatable :: vec_x(:)
            real(dp), allocatable :: vec_y(:)

            call ELL%malloc(4,5,10)
            ELL%data(1,1:3)   = real([9,-3,0],kind=wp)
            ELL%data(2,1:3)   = real([4,7,0],kind=wp)
            ELL%data(3,1:3)   = real([8,-1,8],kind=wp)
            ELL%data(4,1:3)   = real([4,5,6],kind=wp)

            ELL%index(1,1:3) = [1,5,0]
            ELL%index(2,1:3) = [1,2,0]
            ELL%index(3,1:3) = [2,3,4]
            ELL%index(4,1:3) = [1,3,4]
            
            allocate( vec_x(5) , source = 1._wp )
            allocate( vec_y(4) , source = 0._wp )
            call matvec( ELL, vec_x, vec_y )
            
            call check(error, all(vec_y == real([6,11,15,15],kind=wp)) )
            if (allocated(error)) return
        end block
        
    end subroutine

    subroutine test_symmetries(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        block
            integer, parameter :: wp = sp
            type(COO_sp) :: COO
            type(CSR_sp) :: CSR
            real(sp), allocatable :: dense(:,:)
            real(sp), allocatable :: vec_x(:)
            real(sp), allocatable :: vec_y1(:), vec_y2(:), vec_y3(:)

            allocate( vec_x(4)  , source = 1._wp )
            allocate( vec_y1(4) , source = 0._wp )
            allocate( vec_y2(4) , source = 0._wp )
            allocate( vec_y3(4) , source = 0._wp )

            allocate( dense(4,4) , source = &
                    reshape(real([1,0,0,0, &
                                  2,1,0,0, &
                                  0,2,1,0,&
                                  0,0,2,1],kind=wp),[4,4]) )

            call dense2coo( dense , COO )
            COO%sym = k_SYMTRISUP
            call coo2csr(COO, CSR)

            dense(2,1) = 2._wp; dense(3,2) = 2._wp; dense(4,3) = 2._wp
            vec_y1 = matmul( dense, vec_x )
            call check(error, all(vec_y1 == [3,5,5,3]) )
            if (allocated(error)) return

            call matvec( COO , vec_x, vec_y2 )
            call check(error, all(vec_y1 == vec_y2) )
            if (allocated(error)) return

            call matvec( CSR , vec_x, vec_y3 )
            call check(error, all(vec_y1 == vec_y3) )
            if (allocated(error)) return
        end block
        block
            integer, parameter :: wp = dp
            type(COO_dp) :: COO
            type(CSR_dp) :: CSR
            real(dp), allocatable :: dense(:,:)
            real(dp), allocatable :: vec_x(:)
            real(dp), allocatable :: vec_y1(:), vec_y2(:), vec_y3(:)

            allocate( vec_x(4)  , source = 1._wp )
            allocate( vec_y1(4) , source = 0._wp )
            allocate( vec_y2(4) , source = 0._wp )
            allocate( vec_y3(4) , source = 0._wp )

            allocate( dense(4,4) , source = &
                    reshape(real([1,0,0,0, &
                                  2,1,0,0, &
                                  0,2,1,0,&
                                  0,0,2,1],kind=wp),[4,4]) )

            call dense2coo( dense , COO )
            COO%sym = k_SYMTRISUP
            call coo2csr(COO, CSR)

            dense(2,1) = 2._wp; dense(3,2) = 2._wp; dense(4,3) = 2._wp
            vec_y1 = matmul( dense, vec_x )
            call check(error, all(vec_y1 == [3,5,5,3]) )
            if (allocated(error)) return

            call matvec( COO , vec_x, vec_y2 )
            call check(error, all(vec_y1 == vec_y2) )
            if (allocated(error)) return

            call matvec( CSR , vec_x, vec_y3 )
            call check(error, all(vec_y1 == vec_y3) )
            if (allocated(error)) return
        end block
    end subroutine

    subroutine test_cells2sparse(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(COO_dp) :: COO
        integer :: cells(4,2)

        cells(1:4,1) = [2,5,3,4]
        cells(1:4,2) = [4,1,3,2]
        call coo_from_cells(COO,cells,num_points=5,num_cells=2,selfloop=.true.,symtype=k_SYMTRISUP)
        
        call check(error, size(COO%data) == 14 .and. size( COO%index , dim=2 ) == 14 )
        if (allocated(error)) return

    end subroutine

end module test_fsparse