!>---------------------------------------------------
!> Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!>
!> Use of this source code is governed by a MIT
!> license that can be found in the LICENSE.md file
!>---------------------------------------------------
module fsparse_krylov
    use fsparse_constants
    use fsparse_matrix_gallery
    use fsparse_matvec
    use fsparse_conversions
    implicit none
    private

    type, abstract, public :: abs_op_sp
    contains
        procedure(abs_mop_sp), deferred :: times !> matrix operator
        procedure(abs_mop_sp), deferred :: precd !> preconditionner
        procedure(abs_dot_sp), nopass, deferred :: inner !> dot_product
    end type
    type, abstract, public :: abs_op_dp
    contains
        procedure(abs_mop_dp), deferred :: times !> matrix operator
        procedure(abs_mop_dp), deferred :: precd !> preconditionner
        procedure(abs_dot_dp), nopass, deferred :: inner !> dot_product
    end type
    
    abstract interface
        subroutine abs_mop_sp(op,x,y)
            import :: sp, abs_op_sp
            class(abs_op_sp) :: op
            real(sp), intent(in)  :: x(:)
            real(sp), intent(out) :: y(:)
        end subroutine
        real(sp) function abs_dot_sp(x,y) result(r)
            import :: sp
            real(sp), intent(in) :: x(:)
            real(sp), intent(in) :: y(:)
        end function
        subroutine abs_mop_dp(op,x,y)
            import :: dp, abs_op_dp
            class(abs_op_dp) :: op
            real(dp), intent(in)  :: x(:)
            real(dp), intent(out) :: y(:)
        end subroutine
        real(dp) function abs_dot_dp(x,y) result(r)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: y(:)
        end function
    end interface
    
    public :: cgsolve

    interface cgsolve
        module subroutine cgsolve_COO_sp(A,b,x,di,rtol,maxiter,restart)
            type(COO_sp), intent(in), target :: A
            real(sp), intent(in) :: b(:)
            real(sp), intent(out) :: x(:)
            logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
            real(sp), intent(in), optional :: rtol
            integer, intent(in), optional :: maxiter
            logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution
        end subroutine
        module subroutine cgsolve_COO_dp(A,b,x,di,rtol,maxiter,restart)
            type(COO_dp), intent(in), target :: A
            real(dp), intent(in) :: b(:)
            real(dp), intent(out) :: x(:)
            logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
            real(dp), intent(in), optional :: rtol
            integer, intent(in), optional :: maxiter
            logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution
        end subroutine
        module subroutine cgsolve_CSR_sp(A,b,x,di,rtol,maxiter,restart)
            type(CSR_sp), intent(in), target :: A
            real(sp), intent(in) :: b(:)
            real(sp), intent(out) :: x(:)
            logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
            real(sp), intent(in), optional :: rtol
            integer, intent(in), optional :: maxiter
            logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution
        end subroutine
        module subroutine cgsolve_CSR_dp(A,b,x,di,rtol,maxiter,restart)
            type(CSR_dp), intent(in), target :: A
            real(dp), intent(in) :: b(:)
            real(dp), intent(out) :: x(:)
            logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
            real(dp), intent(in), optional :: rtol
            integer, intent(in), optional :: maxiter
            logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution
        end subroutine
        module subroutine cgsolve_dense_sp(A,b,x,di,rtol,maxiter,restart)
            real(sp), intent(in), target :: A(:,:)
            real(sp), intent(in) :: b(:)
            real(sp), intent(out) :: x(:)
            logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
            real(sp), intent(in), optional :: rtol
            integer, intent(in), optional :: maxiter
            logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution
        end subroutine
        module subroutine cgsolve_dense_dp(A,b,x,di,rtol,maxiter,restart)
            real(dp), intent(in), target :: A(:,:)
            real(dp), intent(in) :: b(:)
            real(dp), intent(out) :: x(:)
            logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
            real(dp), intent(in), optional :: rtol
            integer, intent(in), optional :: maxiter
            logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution
        end subroutine
    end interface

    public :: cgsolve_generic
    interface cgsolve_generic
        module procedure :: cgsolve_generic_sp
        module procedure :: cgsolve_generic_dp
    end interface

    contains

    subroutine cgsolve_generic_sp(A,b,x,di,rtol,maxiter,restart)
        class(abs_op_sp), intent(in) :: A
        real(sp), intent(in) :: b(:)
        real(sp), intent(out) :: x(:)
        logical, intent(in) :: di(:)
        integer, intent(in) :: maxiter !! maximum number of iterations
        real(sp), intent(in) :: rtol    !! relative tolerance
        logical, intent(in) :: restart !! use values in x as initial solution

        real(sp), parameter :: zero = zero_sp
        real(sp), parameter :: one  = 1._sp
        real(sp) :: norm_sq, norm_sq0, norm_sq_old, residual0
        real(sp) :: zr1, zr2, zv2, alpha, beta, tolsq
        real(sp), allocatable :: r(:), s(:), p(:), q(:)
        integer :: iter
        !--------------------------------------------
        include 'fsparse_krylov_cgs.inc'
    end subroutine
    subroutine cgsolve_generic_dp(A,b,x,di,rtol,maxiter,restart)
        class(abs_op_dp), intent(in) :: A
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(in) :: di(:)
        integer, intent(in) :: maxiter !! maximum number of iterations
        real(dp), intent(in) :: rtol    !! relative tolerance
        logical, intent(in) :: restart !! use values in x as initial solution

        real(dp), parameter :: zero = zero_dp
        real(dp), parameter :: one  = 1._dp
        real(dp) :: norm_sq, norm_sq0, norm_sq_old, residual0
        real(dp) :: zr1, zr2, zv2, alpha, beta, tolsq
        real(dp), allocatable :: r(:), s(:), p(:), q(:)
        integer :: iter
        !--------------------------------------------
        include 'fsparse_krylov_cgs.inc'
    end subroutine

end module

!======================= submodules for standard types =====================
submodule (fsparse_krylov) fsparse_krylov_CSR
    implicit none
    type, extends(abs_op_sp) :: linop_sp
        type(CSR_sp), pointer :: ptr_CSR_sp ! pointer to the input original matrix
        type(CSR_sp)  :: ptr_pco_CSR_sp ! pointer/memory for out-of-diagonal preconditioner if any
        type(diagonal_sp) :: ptr_pcd_sp    ! pointer/memory for diagonal preconditioner
    contains
        procedure :: times => loc_matvec_sp
        procedure :: precd => loc_precond_sp
        procedure, nopass :: inner => dot_sp
        procedure :: clear => clear_sp
    end type

    type, extends(abs_op_dp) :: linop_dp
        type(CSR_dp), pointer :: ptr_CSR_dp ! pointer to the input original matrix
        type(CSR_dp)  :: ptr_pco_CSR_dp ! pointer/memory for out-of-diagonal preconditioner if any
        type(diagonal_dp) :: ptr_pcd_dp    ! pointer/memory for diagonal preconditioner
    contains
        procedure :: times => loc_matvec_dp
        procedure :: precd => loc_precond_dp
        procedure, nopass :: inner => dot_dp
        procedure :: clear => clear_dp
    end type

    contains
    subroutine loc_matvec_sp(op,x,y)
        class(linop_sp) :: op
        real(sp), intent(in)  :: x(:)
        real(sp), intent(out) :: y(:)
        y = zero_sp
        call matvec(op%ptr_CSR_sp,x,y)
    end subroutine
    subroutine loc_precond_sp(op,x,y)
        class(linop_sp) :: op
        real(sp), intent(in)  :: x(:)
        real(sp), intent(out) :: y(:)
        y = zero_sp
        call matvec(op%ptr_pcd_sp,x,y)
    end subroutine
    real(sp) function dot_sp(x,y) result(r)
        real(sp), intent(in) :: x(:)
        real(sp), intent(in) :: y(:)
        r = dot_product(x,y)
    end function

    subroutine factorization_jacobi_sp(op)
        class(linop_sp) :: op
        allocate(op%ptr_pcd_sp%data(op%ptr_CSR_sp%nrows))
        associate( diag => op%ptr_pcd_sp%data , &
            &       mat => op%ptr_CSR_sp )

            call CSR2diagonal(mat,diag)
            where(abs(diag)>zero_sp) diag = 1._sp/diag

        end associate
    end subroutine

    subroutine clear_sp(op)
        class(linop_sp) :: op
        deallocate(op%ptr_pcd_sp%data)
    end subroutine

    subroutine loc_matvec_dp(op,x,y)
        class(linop_dp) :: op
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y = zero_dp
        call matvec(op%ptr_CSR_dp,x,y)
    end subroutine
    subroutine loc_precond_dp(op,x,y)
        class(linop_dp) :: op
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y = zero_dp
        call matvec(op%ptr_pcd_dp,x,y)
    end subroutine
    real(dp) function dot_dp(x,y) result(r)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)
        r = dot_product(x,y)
    end function

    subroutine factorization_jacobi_dp(op)
        class(linop_dp) :: op
        allocate(op%ptr_pcd_dp%data(op%ptr_CSR_dp%nrows))
        associate( diag => op%ptr_pcd_dp%data , &
            &       mat => op%ptr_CSR_dp )

            call CSR2diagonal(mat,diag)
            where(abs(diag)>zero_dp) diag = 1._dp/diag

        end associate
    end subroutine

    subroutine clear_dp(op)
        class(linop_dp) :: op
        deallocate(op%ptr_pcd_dp%data)
    end subroutine


    module subroutine cgsolve_CSR_sp(A,b,x,di,rtol,maxiter,restart)
        type(CSR_sp), intent(in), target :: A
        real(sp), intent(in) :: b(:) !! right hand side of the linear system
        real(sp), intent(out) :: x(:) !! Solution vector, it can contain also the initial guess
        logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
        real(sp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximun number of iterations
        logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution

        logical, pointer :: di_(:)
        real(sp) :: rtol_
        integer :: maxiter_
        logical :: restart_
        !-------------------------
        type(linop_sp) :: op
        !-------------------------
        rtol_    = 1e-8_sp
        maxiter_ = 1000
        restart_ = .true.
        if(present(maxiter)) maxiter_ = maxiter
        if(present(rtol))    rtol_    = rtol
        if(present(restart)) restart_ = restart
        if(present(di))then
            di_ => di
        else 
            allocate(di_(size(x)),source=.false.)
        end if
        
        !-----------------------------------------------------
        ! set matrix pointer
        op%ptr_CSR_sp => A

        !-----------------------------------------------------
        ! set preconditioner
        call factorization_jacobi_sp(op)

        !-----------------------------------------------------
        ! solve with pccg
        call cgsolve_generic(op,b,x,di_,rtol_,maxiter_,restart_)

        !-----------------------------------------------------
        ! clean memory
        if(.not.present(di)) deallocate(di_)
        call op%clear
    end subroutine

    module subroutine cgsolve_CSR_dp(A,b,x,di,rtol,maxiter,restart)
        type(CSR_dp), intent(in), target :: A
        real(dp), intent(in) :: b(:) !! right hand side of the linear system
        real(dp), intent(out) :: x(:) !! Solution vector, it can contain also the initial guess
        logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximun number of iterations
        logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution

        logical, pointer :: di_(:)
        real(dp) :: rtol_
        integer :: maxiter_
        logical :: restart_
        !-------------------------
        type(linop_dp) :: op
        !-------------------------
        rtol_    = 1e-8_dp
        maxiter_ = 1000
        restart_ = .true.
        if(present(maxiter)) maxiter_ = maxiter
        if(present(rtol))    rtol_    = rtol
        if(present(restart)) restart_ = restart
        if(present(di))then
            di_ => di
        else 
            allocate(di_(size(x)),source=.false.)
        end if
        
        !-----------------------------------------------------
        ! set matrix pointer
        op%ptr_CSR_dp => A

        !-----------------------------------------------------
        ! set preconditioner
        call factorization_jacobi_dp(op)

        !-----------------------------------------------------
        ! solve with pccg
        call cgsolve_generic(op,b,x,di_,rtol_,maxiter_,restart_)

        !-----------------------------------------------------
        ! clean memory
        if(.not.present(di)) deallocate(di_)
        call op%clear
    end subroutine


end submodule
submodule (fsparse_krylov) fsparse_krylov_dense
    implicit none
    type, extends(abs_op_sp) :: linop_sp
        type(dense_sp)  :: ptr_dense_sp ! pointer to the input original matrix
        type(dense_sp)  :: ptr_pco_dense_sp ! pointer/memory for out-of-diagonal preconditioner if any
        type(diagonal_sp) :: ptr_pcd_sp    ! pointer/memory for diagonal preconditioner
    contains
        procedure :: times => loc_matvec_sp
        procedure :: precd => loc_precond_sp
        procedure, nopass :: inner => dot_sp
        procedure :: clear => clear_sp
    end type

    type, extends(abs_op_dp) :: linop_dp
        type(dense_dp)  :: ptr_dense_dp ! pointer to the input original matrix
        type(dense_dp)  :: ptr_pco_dense_dp ! pointer/memory for out-of-diagonal preconditioner if any
        type(diagonal_dp) :: ptr_pcd_dp    ! pointer/memory for diagonal preconditioner
    contains
        procedure :: times => loc_matvec_dp
        procedure :: precd => loc_precond_dp
        procedure, nopass :: inner => dot_dp
        procedure :: clear => clear_dp
    end type

    contains
    subroutine loc_matvec_sp(op,x,y)
        class(linop_sp) :: op
        real(sp), intent(in)  :: x(:)
        real(sp), intent(out) :: y(:)
        y = zero_sp
        call matvec(op%ptr_dense_sp,x,y)
    end subroutine
    subroutine loc_precond_sp(op,x,y)
        class(linop_sp) :: op
        real(sp), intent(in)  :: x(:)
        real(sp), intent(out) :: y(:)
        y = zero_sp
        call matvec(op%ptr_pcd_sp,x,y)
    end subroutine
    real(sp) function dot_sp(x,y) result(r)
        real(sp), intent(in) :: x(:)
        real(sp), intent(in) :: y(:)
        r = dot_product(x,y)
    end function

    subroutine factorization_jacobi_sp(op)
        class(linop_sp) :: op
        allocate(op%ptr_pcd_sp%data(op%ptr_dense_sp%nrows))
        associate( diag => op%ptr_pcd_sp%data , &
            &       mat => op%ptr_dense_sp%data )

            call dense2diagonal(mat,diag)
            where(abs(diag)>zero_sp) diag = 1._sp/diag

        end associate
    end subroutine

    subroutine clear_sp(op)
        class(linop_sp) :: op
        deallocate(op%ptr_pcd_sp%data)
    end subroutine

    subroutine loc_matvec_dp(op,x,y)
        class(linop_dp) :: op
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y = zero_dp
        call matvec(op%ptr_dense_dp,x,y)
    end subroutine
    subroutine loc_precond_dp(op,x,y)
        class(linop_dp) :: op
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y = zero_dp
        call matvec(op%ptr_pcd_dp,x,y)
    end subroutine
    real(dp) function dot_dp(x,y) result(r)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)
        r = dot_product(x,y)
    end function

    subroutine factorization_jacobi_dp(op)
        class(linop_dp) :: op
        allocate(op%ptr_pcd_dp%data(op%ptr_dense_dp%nrows))
        associate( diag => op%ptr_pcd_dp%data , &
            &       mat => op%ptr_dense_dp%data )

            call dense2diagonal(mat,diag)
            where(abs(diag)>zero_dp) diag = 1._dp/diag

        end associate
    end subroutine

    subroutine clear_dp(op)
        class(linop_dp) :: op
        deallocate(op%ptr_pcd_dp%data)
    end subroutine


    module subroutine cgsolve_dense_sp(A,b,x,di,rtol,maxiter,restart)
        real(sp), intent(in), target :: A(:,:)
        real(sp), intent(in) :: b(:) !! right hand side of the linear system
        real(sp), intent(out) :: x(:) !! Solution vector, it can contain also the initial guess
        logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
        real(sp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximun number of iterations
        logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution

        logical, pointer :: di_(:)
        real(sp) :: rtol_
        integer :: maxiter_
        logical :: restart_
        !-------------------------
        type(linop_sp) :: op
        !-------------------------
        rtol_    = 1e-8_sp
        maxiter_ = 1000
        restart_ = .true.
        if(present(maxiter)) maxiter_ = maxiter
        if(present(rtol))    rtol_    = rtol
        if(present(restart)) restart_ = restart
        if(present(di))then
            di_ => di
        else 
            allocate(di_(size(x)),source=.false.)
        end if
        
        !-----------------------------------------------------
        ! set matrix pointer
        op%ptr_dense_sp%data => A(:,:)

        !-----------------------------------------------------
        ! set preconditioner
        call factorization_jacobi_sp(op)

        !-----------------------------------------------------
        ! solve with pccg
        call cgsolve_generic(op,b,x,di_,rtol_,maxiter_,restart_)

        !-----------------------------------------------------
        ! clean memory
        if(.not.present(di)) deallocate(di_)
        call op%clear
    end subroutine

    module subroutine cgsolve_dense_dp(A,b,x,di,rtol,maxiter,restart)
        real(dp), intent(in), target :: A(:,:)
        real(dp), intent(in) :: b(:) !! right hand side of the linear system
        real(dp), intent(out) :: x(:) !! Solution vector, it can contain also the initial guess
        logical, intent(in), optional, target :: di(:) !! Dirichlet conditions for blocked d.o.f
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximun number of iterations
        logical, intent(in), optional :: restart !! Set X=0 or use it as starting solution

        logical, pointer :: di_(:)
        real(dp) :: rtol_
        integer :: maxiter_
        logical :: restart_
        !-------------------------
        type(linop_dp) :: op
        !-------------------------
        rtol_    = 1e-8_dp
        maxiter_ = 1000
        restart_ = .true.
        if(present(maxiter)) maxiter_ = maxiter
        if(present(rtol))    rtol_    = rtol
        if(present(restart)) restart_ = restart
        if(present(di))then
            di_ => di
        else 
            allocate(di_(size(x)),source=.false.)
        end if
        
        !-----------------------------------------------------
        ! set matrix pointer
        op%ptr_dense_dp%data => A(:,:)

        !-----------------------------------------------------
        ! set preconditioner
        call factorization_jacobi_dp(op)

        !-----------------------------------------------------
        ! solve with pccg
        call cgsolve_generic(op,b,x,di_,rtol_,maxiter_,restart_)

        !-----------------------------------------------------
        ! clean memory
        if(.not.present(di)) deallocate(di_)
        call op%clear
    end subroutine


end submodule

