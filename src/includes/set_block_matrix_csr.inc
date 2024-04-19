        integer :: i, j, k, ik, jk, adr1, adr2
        
        select case(self%sym)
        case(k_NOSYMMETRY)
            do i = 1 , num_dof
                ik = adr_dof(i)
                do j = 1 , num_dof
                    jk = adr_dof(j)
                    do k = self%rowptr(ik), self%rowptr(ik+1)-1
                        if( jk == self%col(k) ) then
                            self%data(k) = self%data(k) + block_matrix(i,j)
                            exit
                        end if
                    end do
                end do
            end do
        case(k_SYMTRIINF)
			do i = 1 , num_dof
				ik = adr_dof(i)
				k = self%rowptr(ik+1)-1
				self%data(k) = self%data(k) + block_matrix(i,i)
                do j = 1 , i-1
                    jk = adr_dof(j)
					adr1 = max(ik,jk)
					adr2 = min(ik,jk)
					k = self%rowptr(adr1)
                    do while ( self%col(k) .ne. adr2 )
						k = k + 1
                    end do
					self%data(k) = self%data(k) + block_matrix(j,i)
                end do
            end do

        case(k_SYMTRISUP)
            do i = 1 , num_dof
				ik = adr_dof(i)
				k = self%rowptr(ik)
				self%data(k) = self%data(k) + block_matrix(i,i)
                do j = i+1 , num_dof
                    jk = adr_dof(j)
					adr1 = max(ik,jk)
					adr2 = min(ik,jk)
					k = self%rowptr(adr1)
                    do while ( self%col(k) .ne. adr2 )
						k = k + 1
                    end do
					self%data(k) = self%data(k) + block_matrix(j,i)
                end do
            end do

        end select