!---------------------------------------------------
! Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------

module fsparse_cells2sparse
    use fsparse_matrix_gallery
    use fsparse_sort, only: coo2ordered
    implicit none
    private
    public :: coo_from_cells
    
contains
    
    subroutine coo_from_cells(COO,cells,num_points,num_cells,shift,symtype,selfloop)
        class(COO_t), intent(inout) :: COO
        integer, intent(in) :: cells(:,:) !! Cell to point connectivity (a cell can be a volume element, a face, an edge, etc)
        integer, intent(in) :: num_points !! total number of points grouped by the cells
        integer, intent(in), optional :: num_cells !! Number of cells, if not given then computed from the size of the cells array
        integer, intent(in), optional :: shift    !! Shift the numbering of points
        integer, intent(in), optional :: symtype  !! 0=> full matrix, 1=> triangular inf, 2=> triangular sup
        logical, intent(in), optional :: selfloop !! include diagonal elements
        ! -- Internal Variables
        integer :: i, j, n_i, n_j, adr, adr0, adrg, ned, num_edges_max, cell
        integer :: npc !! number of points per cell
        integer :: num_cells_, shift_, symtype_
        logical :: selfloop_
        !---------------------------------------------------------
        num_cells_ = 0
        shift_     = 0
        symtype_   = k_NOSYMMETRY
        selfloop_  = .true.
        if(present(num_cells)) then
            num_cells_  = num_cells
        else
            num_cells_  = size( cells , dim=2 )
        end if
        if(present(shift))    shift_    = shift
        if(present(symtype))  symtype_  = symtype
        if(present(selfloop)) selfloop_ = selfloop
        COO%sym = symtype_
        npc = size( cells , dim=1 )
        adr0 = COO%nnz
        !---------------------------------------------------------
        ! determine array sizes for pre-allocation
        ned = int( npc*(npc-1)/2 ) !! number of edges per element
        if (symtype_ == k_NOSYMMETRY) ned = 2*ned
            
        num_edges_max = int(npc*(npc-1)/2)*num_cells_
        if (symtype_ == k_NOSYMMETRY) num_edges_max = 2*num_edges_max
        if (selfloop_) num_edges_max = num_edges_max + num_points
        !---------------------------------------------------------
        ! (re)allocate memory
        call coo%malloc( num_points+shift_ , num_points+shift_ , num_edges_max )
        !---------------------------------------------------------
        ! Fill in COO sparse data
        if( selfloop_ ) then ! diagonal terms
        COO%index(1,adr0+1:adr0+num_points) = [ (i, i = shift_+1, shift_+num_points) ]
        COO%index(2,adr0+1:adr0+num_points) = COO%index(1,adr0+1:adr0+num_points)
        adr0 = adr0 + num_points
        end if
            
        select case( symtype_ ) ! off-diagonal terms
        case( k_NOSYMMETRY ) ! Full matrix
            do concurrent (cell = 1:num_cells_)
                adrg = adr0 + ned*(cell-1)
                do i = 1, npc
                    do j = 1, i-1
                        adr = 2*int((i-1)*(i-2)/2)+2*(j-1)+1
                        n_i = cells(i,cell)+shift_
                        n_j = cells(j,cell)+shift_
                        COO%index(1:2,adrg+adr)   = [ n_i , n_j ]
                        COO%index(1:2,adrg+adr+1) = [ n_j , n_i ]
                    end do
                end do
            end do
            
        case( k_SYMTRIINF ) ! Triangular inferior storage
            do concurrent (cell = 1:num_cells_)
                adrg = adr0 + ned*(cell-1)
                do i = 1, npc
                    do j = 1, i-1
                        adr = int((i-1)*(i-2)/2) + j
                        n_i = cells(i,cell)+shift_
                        n_j = cells(j,cell)+shift_
                        COO%index(1:2,adrg+adr) = [ max(n_i,n_j) , min(n_i,n_j) ]
                    end do
                end do
            end do
            
        case( k_SYMTRISUP ) ! Triangular superior storage
            do concurrent (cell = 1:num_cells_)
                adrg = adr0 + ned*(cell-1)
                do i = 1, npc
                    do j = 1, i-1
                        adr = int((i-1)*(i-2)/2) + j
                        n_i = cells(i,cell)+shift_
                        n_j = cells(j,cell)+shift_
                        COO%index(1:2,adrg+adr) = [ min(n_i,n_j) , max(n_i,n_j) ]
                    end do
                end do
            end do
            
        end select
        
        COO%isOrdered = .false.
        call coo2ordered(COO)
    end subroutine
    
end module fsparse_cells2sparse