function x = get_solution(A, m,n, indx)
% this routine provided the matrix M(m x m+1) obtained
% by calling ReducedRowEchelonForm(..) and corresponding indx vector
% returns to x the solution of the system Ax=b.
    for i=1:m
        x(indx(i),1) = A(i,n);
    end
end

