!---------------------------------------------------
! Copyright 2023-present Transvalor S.A.
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module spmv_kernels
    use iso_fortran_env
    implicit none
    private

    interface spmv_coo
        module procedure spmv_coo_1d_sp
        module procedure spmv_coo_1d_dp
        module procedure spmv_coo_2d_sp
        module procedure spmv_coo_2d_dp
    end interface

    interface spmv_coo_sym
        module procedure spmv_coo_sym_1d_sp
        module procedure spmv_coo_sym_1d_dp
        module procedure spmv_coo_sym_2d_sp
        module procedure spmv_coo_sym_2d_dp
    end interface

    interface spmv_csr
        module procedure spmv_csr_1d_sp
        module procedure spmv_csr_1d_dp
        module procedure spmv_csr_2d_sp
        module procedure spmv_csr_2d_dp
    end interface

    interface spmv_csc
        module procedure spmv_csc_1d_sp
        module procedure spmv_csc_1d_dp
        module procedure spmv_csc_2d_sp
        module procedure spmv_csc_2d_dp
    end interface

    interface spmv_ell
        module procedure spmv_ell_1d_sp
        module procedure spmv_ell_1d_dp
        module procedure spmv_ell_2d_sp
        module procedure spmv_ell_2d_dp
    end interface
        
    public :: spmv_coo, spmv_coo_sym
    public :: spmv_csr
    public :: spmv_csc
    public :: spmv_ell
    contains

    !! spmv_coo_kernels
    subroutine spmv_coo_1d_sp(data,index,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_coo_1d.inc'
    end subroutine

    subroutine spmv_coo_1d_dp(data,index,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_coo_1d.inc'
    end subroutine

    subroutine spmv_coo_2d_sp(data,index,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_coo_2d.inc'
    end subroutine

    subroutine spmv_coo_2d_dp(data,index,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_coo_2d.inc'
    end subroutine

    subroutine spmv_coo_sym_1d_sp(data,index,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_coo_sym_1d.inc'
    end subroutine

    subroutine spmv_coo_sym_1d_dp(data,index,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_coo_sym_1d.inc'
    end subroutine

    subroutine spmv_coo_sym_2d_sp(data,index,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_coo_sym_2d.inc'
    end subroutine

    subroutine spmv_coo_sym_2d_dp(data,index,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_coo_sym_2d.inc'
    end subroutine

    !! spmv_csr_kernels
    subroutine spmv_csr_1d_sp(data,col,rowptr,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_csr_1d.inc'
    end subroutine

    subroutine spmv_csr_1d_dp(data,col,rowptr,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_csr_1d.inc'
    end subroutine

    subroutine spmv_csr_2d_sp(data,col,rowptr,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_csr_2d.inc'
    end subroutine

    subroutine spmv_csr_2d_dp(data,col,rowptr,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_csr_2d.inc'
    end subroutine

    !! spmv_csc_kernels
    subroutine spmv_csc_1d_sp(data,row,colptr,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_csc_1d.inc'
    end subroutine

    subroutine spmv_csc_1d_dp(data,row,colptr,N,M,NNZ,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_csc_1d.inc'
    end subroutine

    subroutine spmv_csc_2d_sp(data,row,colptr,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_csc_2d.inc'
    end subroutine

    subroutine spmv_csc_2d_dp(data,row,colptr,N,M,NNZ,dim,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_csc_2d.inc'
    end subroutine

    !! spmv_ell_kernels
    subroutine spmv_ell_1d_sp(data,index,N,M,MNZ_P_ROW,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_ell_1d.inc'
    end subroutine

    subroutine spmv_ell_1d_dp(data,index,N,M,MNZ_P_ROW,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_ell_1d.inc'
    end subroutine

    subroutine spmv_ell_2d_sp(data,index,N,M,MNZ_P_ROW,dim,vec_x,vec_y)
        integer, parameter      :: wp = real32
        include 'spmv_ell_2d.inc'
    end subroutine

    subroutine spmv_ell_2d_dp(data,index,N,M,MNZ_P_ROW,dim,vec_x,vec_y)
        integer, parameter      :: wp = real64
        include 'spmv_ell_2d.inc'
    end subroutine

end module