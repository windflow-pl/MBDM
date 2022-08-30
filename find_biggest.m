function b = find_biggest(M,m,n,k,s)
    big = 0;
    b=[0,0];
    if (s==1)
        for j=k:n-1
            for i=k:m
                if(big<abs(M(i,j)))
                    big = abs(M(i,j));
                    b(1)=i;
                    b(2)=j;
                end
            end
        end
    end
    if (s==0)
        for j=k:n
            for i=k:m
                if(big<abs(M(i,j)))
                    big = abs(M(i,j));
                    b(1)=i;
                    b(2)=j;
                end
            end
        end
    end 
end