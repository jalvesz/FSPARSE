!---------------------------------------------------
! Copyright 2023-present Transvalor S.A.
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module fsparse_matvec
    use fsparse_matrix_gallery
    use spmv_kernels
    use iso_fortran_env
    implicit none
    private

    interface matvec
        module procedure matvec_coo_1d_sp
        module procedure matvec_coo_1d_dp
        module procedure matvec_coo_2d_sp
        module procedure matvec_coo_2d_dp

        module procedure matvec_csr_1d_sp
        module procedure matvec_csr_1d_dp
        module procedure matvec_csr_2d_sp
        module procedure matvec_csr_2d_dp

        module procedure matvec_csc_1d_sp
        module procedure matvec_csc_1d_dp
        module procedure matvec_csc_2d_sp
        module procedure matvec_csc_2d_dp

        module procedure matvec_ell_1d_sp
        module procedure matvec_ell_1d_dp
        module procedure matvec_ell_2d_sp
        module procedure matvec_ell_2d_dp
    end interface

    public :: matvec
    contains
    
    !! matvec_coo
    subroutine matvec_coo_1d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(COOr32_t), intent(in) :: matrix
        include 'matvec_coo_1d.inc'
    end subroutine

    subroutine matvec_coo_1d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(COOr64_t), intent(in) :: matrix
        include 'matvec_coo_1d.inc'
    end subroutine

    subroutine matvec_coo_2d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(COOr32_t), intent(in) :: matrix
        include 'matvec_coo_2d.inc'
    end subroutine

    subroutine matvec_coo_2d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(COOr64_t), intent(in) :: matrix
        include 'matvec_coo_2d.inc'
    end subroutine

    !! matvec_csr
    subroutine matvec_csr_1d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSRr32_t), intent(in) :: matrix
        include 'matvec_csr_1d.inc'
    end subroutine

    subroutine matvec_csr_1d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSRr64_t), intent(in) :: matrix
        include 'matvec_csr_1d.inc'
    end subroutine

    subroutine matvec_csr_2d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSRr32_t), intent(in) :: matrix
        include 'matvec_csr_2d.inc'
    end subroutine

    subroutine matvec_csr_2d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSRr64_t), intent(in) :: matrix
        include 'matvec_csr_2d.inc'
    end subroutine

    !! matvec_csc
    subroutine matvec_csc_1d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSCr32_t), intent(in) :: matrix
        include 'matvec_csc_1d.inc'
    end subroutine

    subroutine matvec_csc_1d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSCr64_t), intent(in) :: matrix
        include 'matvec_csc_1d.inc'
    end subroutine

    subroutine matvec_csc_2d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSCr32_t), intent(in) :: matrix
        include 'matvec_csc_2d.inc'
    end subroutine

    subroutine matvec_csc_2d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSCr64_t), intent(in) :: matrix
        include 'matvec_csc_2d.inc'
    end subroutine

    !! matvec_ell
    subroutine matvec_ell_1d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(ELLr32_t), intent(in) :: matrix
        include 'matvec_ell_1d.inc'
    end subroutine

    subroutine matvec_ell_1d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(ELLr64_t), intent(in) :: matrix
        include 'matvec_ell_1d.inc'
    end subroutine

    subroutine matvec_ell_2d_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(ELLr32_t), intent(in) :: matrix
        include 'matvec_ell_2d.inc'
    end subroutine

    subroutine matvec_ell_2d_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(ELLr64_t), intent(in) :: matrix
        include 'matvec_ell_2d.inc'
    end subroutine

end module fsparse_matvec