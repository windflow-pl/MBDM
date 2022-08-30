function  [y_out, er] = Num_Int(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,h)
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
    % store input q and q_prime
    q_n = q;
    q_prime_n = q_prime;
    % =====================================================================
    % ======= MAIN ROUTINE
    % =====================================================================
    % == k1
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,h);
    if er==1
        fprintf('error while calculating k1 for R-K\n');
        y_out=[0;0];
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
    f = [v_p_out;v_pp_out];
    k1 = h*f;
    % == k2
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = q_n(i,1) + 0.5*k1(m+1,1);
            q_prime(i,1) = q_prime_n(i,1) + 0.5*k1(nind + m+1,1);
        end
    end
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t+0.5*h,h);
    if er==1
        fprintf('error while calculating k2 for R-K\n');
        y_out=0;
        return;
    end
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f = [v_p_out;v_pp_out];
    k2 = h*f;
    % == k3
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = q_n(i,1) + 0.5*k2(m+1,1);
            q_prime(i,1) = q_prime_n(i,1) + 0.5*k2(nind + m+1,1);
        end
    end
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t+0.5*h,h);
    if er==1
        fprintf('error while calculating k3 for R-K\n');
        y_out=0;
        return;
    end
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f = [v_p_out;v_pp_out];
    k3 = h*f;
    % == k4
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = q_n(i,1)+k3(m+1,1);
            q_prime(i,1) = q_prime_n(i,1)+k3(nind + m+1,1);
        end
    end
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t+h,h);
    if er==1
        fprintf('error while calculating k4 for R-K\n');
        y_out=0;
        return;
    end
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f = [v_p_out;v_pp_out];
    k4 = h*f;
    % == y_out
    y_out = y+(1/6)*(k1+2*k2+2*k3+k4);
% %     fprintf('================\n');
%       y_out = y + k1;

end
