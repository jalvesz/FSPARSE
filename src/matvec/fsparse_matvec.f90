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
        module procedure matvec_coo_1d_r_sp
        module procedure matvec_coo_1d_r_dp
        module procedure matvec_coo_2d_r_sp
        module procedure matvec_coo_2d_r_dp
        module procedure matvec_coo_1d_c_sp
        module procedure matvec_coo_1d_c_dp
        module procedure matvec_coo_2d_c_sp
        module procedure matvec_coo_2d_c_dp

        module procedure matvec_csr_1d_r_sp
        module procedure matvec_csr_1d_r_dp
        module procedure matvec_csr_2d_r_sp
        module procedure matvec_csr_2d_r_dp
        module procedure matvec_csr_1d_c_sp
        module procedure matvec_csr_1d_c_dp
        module procedure matvec_csr_2d_c_sp
        module procedure matvec_csr_2d_c_dp

        module procedure matvec_csc_1d_r_sp
        module procedure matvec_csc_1d_r_dp
        module procedure matvec_csc_2d_r_sp
        module procedure matvec_csc_2d_r_dp
        module procedure matvec_csc_1d_c_sp
        module procedure matvec_csc_1d_c_dp
        module procedure matvec_csc_2d_c_sp
        module procedure matvec_csc_2d_c_dp

        module procedure matvec_ell_1d_r_sp
        module procedure matvec_ell_1d_r_dp
        module procedure matvec_ell_2d_r_sp
        module procedure matvec_ell_2d_r_dp
        module procedure matvec_ell_1d_c_sp
        module procedure matvec_ell_1d_c_dp
        module procedure matvec_ell_2d_c_sp
        module procedure matvec_ell_2d_c_dp
    end interface

    public :: matvec
    contains
    
    !! matvec_coo
    subroutine matvec_coo_1d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(COOr32_t), intent(in) :: matrix
        include 'matvec_coo_r_1d.inc'
    end subroutine

    subroutine matvec_coo_1d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(COOr64_t), intent(in) :: matrix
        include 'matvec_coo_r_1d.inc'
    end subroutine

    subroutine matvec_coo_2d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(COOr32_t), intent(in) :: matrix
        include 'matvec_coo_r_2d.inc'
    end subroutine

    subroutine matvec_coo_2d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(COOr64_t), intent(in) :: matrix
        include 'matvec_coo_r_2d.inc'
    end subroutine
  
    subroutine matvec_coo_1d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(COOc32_t), intent(in) :: matrix
        include 'matvec_coo_c_1d.inc'
    end subroutine

    subroutine matvec_coo_1d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(COOc64_t), intent(in) :: matrix
        include 'matvec_coo_c_1d.inc'
    end subroutine

    subroutine matvec_coo_2d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(COOc32_t), intent(in) :: matrix
        include 'matvec_coo_c_2d.inc'
    end subroutine

    subroutine matvec_coo_2d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(COOc64_t), intent(in) :: matrix
        include 'matvec_coo_c_2d.inc'
    end subroutine


    
    !! matvec_csr
    subroutine matvec_csr_1d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSRr32_t), intent(in) :: matrix
        include 'matvec_csr_r_1d.inc'
    end subroutine

    subroutine matvec_csr_1d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSRr64_t), intent(in) :: matrix
        include 'matvec_csr_r_1d.inc'
    end subroutine

    subroutine matvec_csr_2d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSRr32_t), intent(in) :: matrix
        include 'matvec_csr_r_2d.inc'
    end subroutine

    subroutine matvec_csr_2d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSRr64_t), intent(in) :: matrix
        include 'matvec_csr_r_2d.inc'
    end subroutine

    subroutine matvec_csr_1d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSRc32_t), intent(in) :: matrix
        include 'matvec_csr_c_1d.inc'
    end subroutine

    subroutine matvec_csr_1d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSRc64_t), intent(in) :: matrix
        include 'matvec_csr_c_1d.inc'
    end subroutine

    subroutine matvec_csr_2d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSRc32_t), intent(in) :: matrix
        include 'matvec_csr_c_2d.inc'
    end subroutine

    subroutine matvec_csr_2d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSRc64_t), intent(in) :: matrix
        include 'matvec_csr_c_2d.inc'
    end subroutine

    !! matvec_csc
    subroutine matvec_csc_1d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSCr32_t), intent(in) :: matrix
        include 'matvec_csc_r_1d.inc'
    end subroutine

    subroutine matvec_csc_1d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSCr64_t), intent(in) :: matrix
        include 'matvec_csc_r_1d.inc'
    end subroutine

    subroutine matvec_csc_2d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSCr32_t), intent(in) :: matrix
        include 'matvec_csc_r_2d.inc'
    end subroutine

    subroutine matvec_csc_2d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSCr64_t), intent(in) :: matrix
        include 'matvec_csc_r_2d.inc'
    end subroutine

    subroutine matvec_csc_1d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSCc32_t), intent(in) :: matrix
        include 'matvec_csc_c_1d.inc'
    end subroutine

    subroutine matvec_csc_1d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSCc64_t), intent(in) :: matrix
        include 'matvec_csc_c_1d.inc'
    end subroutine

    subroutine matvec_csc_2d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(CSCc32_t), intent(in) :: matrix
        include 'matvec_csc_c_2d.inc'
    end subroutine

    subroutine matvec_csc_2d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(CSCc64_t), intent(in) :: matrix
        include 'matvec_csc_c_2d.inc'
    end subroutine

    !! matvec_ell
    subroutine matvec_ell_1d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(ELLr32_t), intent(in) :: matrix
        include 'matvec_ell_r_1d.inc'
    end subroutine

    subroutine matvec_ell_1d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(ELLr64_t), intent(in) :: matrix
        include 'matvec_ell_r_1d.inc'
    end subroutine

    subroutine matvec_ell_2d_r_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(ELLr32_t), intent(in) :: matrix
        include 'matvec_ell_r_2d.inc'
    end subroutine

    subroutine matvec_ell_2d_r_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(ELLr64_t), intent(in) :: matrix
        include 'matvec_ell_r_2d.inc'
    end subroutine

    subroutine matvec_ell_1d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(ELLc32_t), intent(in) :: matrix
        include 'matvec_ell_c_1d.inc'
    end subroutine

    subroutine matvec_ell_1d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(ELLc64_t), intent(in) :: matrix
        include 'matvec_ell_c_1d.inc'
    end subroutine

    subroutine matvec_ell_2d_c_sp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real32
        type(ELLc32_t), intent(in) :: matrix
        include 'matvec_ell_c_2d.inc'
    end subroutine

    subroutine matvec_ell_2d_c_dp(matrix,vec_x,vec_y)
        integer, parameter      :: wp = real64
        type(ELLc64_t), intent(in) :: matrix
        include 'matvec_ell_c_2d.inc'
    end subroutine

end module fsparse_matvec
