!>---------------------------------------------------
!> Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!>
!> Use of this source code is governed by a MIT
!> license that can be found in the LICENSE.md file
!>---------------------------------------------------
module test_solvers
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fsparse_krylov
    implicit none
    private
    public :: collect_suite
    contains

    !> Collect all exported unit tests
    subroutine collect_suite(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest('pccg', test_pccg) &
        ]
    end subroutine

    subroutine test_pccg(error)
        !> Error handling
        type(error_type), allocatable, intent(out) :: error
        real :: A(2,2) = reshape([4.0,1.0,1.0,3.0],shape=[2,2])
        real :: b(2) = [1.0,2.0]
        real :: x(2)
        real :: tol = 0.0001

        call cgsolve(A,b,x)
        
        call check(error, all(abs(x - [1.0/11,7.0/11])<epsilon(0.0)) )
        if (allocated(error)) return
    end subroutine

end module test_solvers