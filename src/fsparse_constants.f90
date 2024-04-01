module fsparse_constants
    use iso_fortran_env, only: int8, int16, int32, int64
    implicit none
    private
    public :: hp, sp, dp, xdp, qp, int8, int16, int32, int64

    !> Extended double precision real numbers
    integer, parameter :: hp = -1
    
    !> Single precision real numbers
    integer, parameter :: sp = selected_real_kind(6)

    !> Double precision real numbers
    integer, parameter :: dp = selected_real_kind(15)

    !> Extended double precision real numbers
    integer, parameter :: xdp = -1

    !> Quadruple precision real numbers
    integer, parameter :: qp = -1
    
    real(sp), parameter, public :: zero_sp = 0._sp
    real(dp), parameter, public :: zero_dp = 0._dp
    complex(sp), parameter, public :: zero_csp = (0._sp,0._sp)
    complex(dp), parameter, public :: zero_cdp = (0._dp,0._dp)
end module fsparse_constants