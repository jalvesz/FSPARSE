!---------------------------------------------------
! Copyright 2023-present Transvalor S.A. (José R. Alves Z.)
!
! Use of this source code is governed by a MIT
! license that can be found in the LICENSE.md file
!---------------------------------------------------
module fsparse
  !! User API: by loading this module all matrix types and methods are directly referenced
  use fsparse_matrix_gallery
  use fsparse_matvec
  use fsparse_sort
  use fsparse_conversions
  use fsparse_cells2sparse
end module fsparse