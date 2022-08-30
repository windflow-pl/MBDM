function [A,indx,er] = RowReducedEchelonForm(M,s)
% returns to M(m x n) its reduced row echelon form
% s indicates the solution, and if s==1 it will turn to solver
% of the system in form of Ax=b; where M is expected to be
% M(m x m+1)
% if s!=1 the routine will follow as usual to get reduced row echelon form of M.
% Routine is based on Gauss ellimination
    A = M;
    [m,n] = size(A);
    er =0;
    for i=1:n
        indx(i) = i;
    end
    for i=1:m
        b = find_biggest(A,m,n,i,s);
        if (b == [0,0])
            disp(sprintf('Singular matrix in routine RowReducedEchelonForm.\n'));
            er = 1;
            break;
        end
        if (A(b(1),b(2))==0)
            disp(sprintf('Singular matrix in routine RowReducedEchelonForm\n'));
            er = 1;
            break;
        end
        A([i b(1)],:) = A([b(1) i],:);% swap rows
        A(:,[i b(2)]) = A(:,[b(2) i]);% swap column
        tmp = indx(i); indx(i) = indx(b(2)); indx(b(2)) = tmp;
        clear tmp;
        div = A(i,i);
        for j=1:n
            A(i,j) = A(i,j)/div;
        end
        for k=1:m
            if(k~=i)
                div = A(k,i);
                for j=i:n
                    A(k,j) = A(k,j) - div*A(i,j);
                end
            end
        end

    end
    
    

%         div = M[i][i];
%         for(j=i;j<n;j++){
%             M[i][j]=M[i][j]/div;
%         }
%         for(k=0;k<m;k++){
%             if(k!=i){
%                 div = M[k][i];
%                 for(j=i;j<n;j++){
%                     M[k][j] = M[k][j] - div*M[i][j];
%                 }
%             }
%         }
%     }
%     free(b);



end


                