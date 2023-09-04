program test_main
use test_fsparse
implicit none

call test_coo()
call test_csr()
call test_csc()
call test_ell()

end program
