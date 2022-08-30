function [FA] = calc_FA(t,q,q_prime,nBody,body)
    % allocate
    FA = zeros(3*nBody,1);
    % ==== add gravity
    for i=1:nBody
        m=(i-1)*3;
        FA(m+1:m+3,1) = [0;0;-body(i).m*9.81];
    end
end
