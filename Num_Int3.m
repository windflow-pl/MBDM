function  [y_out, er] = Num_Int3(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,h)
    % number of independent variables
    nind = 7*nBody-nCon;
    % allocate v and v_p
    v = zeros(7*nBody-nCon,1);
    v_p = zeros(7*nBody-nCon,1);
    % extract v
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v(m+1,1) = q(i,1);
            v_p(m+1,1) = q_prime(i,1);
            m=m+1;
        end
    end
    % assign y
    y = zeros(2*(7*nBody-nCon),1);
    y = [v;v_p];
    y_out=zeros(2*(7*nBody-nCon),1);
    % store input q and q_prime
    q_n = q;
    q_prime_n = q_prime;
    % =====================================================================
    % ======= MAIN ROUTINE (NUMERICAL INTEGRATION ROUTINE)
    % =====================================================================
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,h);
    if er==1
        fprintf('error while calculating k1 for R-K\n');
        return;        
    end
    % extract v_p_out and v_pp_out
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    y_p = [v_p_out;v_pp_out];
    y_out = y + y_p*h;
%     y
%     y_p
%     n
    
    
end