allocate(R,mold=B)
allocate(S,mold=B)
allocate(P,mold=B)
allocate(Q,mold=B)

norm_sq = A%inner( B, B )
norm_sq0 = norm_sq
residual0 = sqrt(norm_sq0)
if ( norm_sq0 > zero ) then
    if(restart) X = zero
    where ( di ) X = B
    
    call A%times(X, R) !< Apply matrix operator
    R = B - R
    where ( di ) R = zero
    
    call A%precd(R,P) !< Precondition Residual
    
    where ( di ) P = zero
    
    tolsq = rtol*rtol
    iter = 0
    zr1 = zero
    zr2 = one
    do while ( (iter < maxiter) .AND. (norm_sq > tolsq * norm_sq0) )
        
        call A%precd(R,S)
        where ( di ) S = zero
        zr2 = A%inner( R, S )
    
        if (iter>0) then
            beta = zr2 / zr1
            P = S + beta * P
        end if
        
        call A%times(P, Q)
        where ( di ) Q = zero
        zv2 = A%inner( P, Q )
    
        alpha = zr2 / zv2
    
        X = X + alpha * P
        R = R - alpha * Q
        norm_sq = A%inner( R, R )
        norm_sq_old = norm_sq
        zr1 = zr2
        iter = iter + 1
    end do
end if