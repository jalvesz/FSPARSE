module test_fsparse
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
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
            new_unittest('symmetries', test_symmetries) &
        ]
    end subroutine

    subroutine test_coo(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error

        real(sp), allocatable :: dense(:,:)
        type(COOr32_t) :: COO
        real(sp), allocatable :: vec_x(:)
        real(sp), allocatable :: vec_y1(:), vec_y2(:)

        allocate( dense(4,5) , source = &
                reshape([9.0,4.0,0.0,4.0, &
                         0.0,7.0,8.0,0.0, &
                         0.0,0.0,-1.0,5.0,&
                         0.0,0.0,8.0,6.0,&
                         -3.0,0.0,0.0,0.0],[4,5]) )
        
        call dense2coo( dense , COO )
        
        allocate( vec_x(5) , source = 1._sp )
        allocate( vec_y1(4) , source = 0._sp )
        allocate( vec_y2(4) , source = 0._sp )

        vec_y1 = matmul( dense, vec_x )
        
        call check(error, all(vec_y1 == [6.0,11.0,15.0,15.0]) )
        if (allocated(error)) return

        call matvec( COO , vec_x, vec_y2 )
        call check(error, all(vec_y1 == vec_y2) )
        if (allocated(error)) return
        
    end subroutine

    subroutine test_coo2ordered(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        type(COOr32_t) :: COO

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
        
        call check(error, all(vec_y == dble([6.0,11.0,15.0,15.0])) )
        if (allocated(error)) return
        
    end subroutine

    subroutine test_csc(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
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

        call check(error, all(vec_y == dble([6.0,11.0,15.0,15.0])) )
        if (allocated(error)) return
        
    end subroutine

    subroutine test_ell(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
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

        call check(error, all(vec_y == [6.0,11.0,15.0,15.0]) )
        if (allocated(error)) return
        
    end subroutine

    subroutine test_symmetries(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error

        real(sp), allocatable :: dense(:,:)
        type(COOr32_t) :: COO
        type(CSRr32_t) :: CSR
        real(sp), allocatable :: vec_x(:)
        real(sp), allocatable :: vec_y1(:), vec_y2(:), vec_y3(:)

        allocate( vec_x(4) , source = 1._sp )
        allocate( vec_y1(4) , source = 0._sp )
        allocate( vec_y2(4) , source = 0._sp )
        allocate( vec_y3(4) , source = 0._sp )

        allocate( dense(4,4) , source = &
                reshape([1.0,0.0,0.0,0.0, &
                         2.0,1.0,0.0,0.0, &
                         0.0,2.0,1.0,0.0,&
                         0.0,0.0,2.0,1.0],[4,4]) )

        call dense2coo( dense , COO )
        COO%sym = k_SYMTRISUP

        dense(2,1) = 2.0; dense(3,2) = 2.0; dense(4,3) = 2.0
        vec_y1 = matmul( dense, vec_x )

        call check(error, all(vec_y1 == [3.0,5.0,5.0,3.0]) )
        if (allocated(error)) return

        call matvec( COO , vec_x, vec_y2 )
        call check(error, all(vec_y1 == vec_y2) )
        if (allocated(error)) return

        call coo2csr(COO, CSR)
        call matvec( CSR , vec_x, vec_y3 )
        call check(error, all(vec_y1 == vec_y3) )
        if (allocated(error)) return

    end subroutine

end module test_fsparse

program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use test_fsparse, only : collect_suite
    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
  
    stat = 0
  
    testsuites = [ &
      new_testsuite("fsparse", collect_suite) &
      ]
  
    do is = 1, size(testsuites)
      write(error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  
    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
    end if
  
  end program tester