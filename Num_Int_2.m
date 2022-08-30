function  [y_out, er] = Num_Int_2(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,h,q_store,q_prime_store)
    % number of independent variables
    nind = 7*nBody-nCon;
    % allocate v and v_p
    v = zeros(7*nBody-nCon,1);
    v_p = zeros(7*nBody-nCon,1);
    % extract v and v_p
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v(m+1,1) = q(i,1);
            v_p(m+1,1) = q_prime(i,1);
            m=m+1;
        end
    end
    % assign y
    y_n = zeros(2*(7*nBody-nCon),1);
    y_n = [v;v_p];
    % store input q and q_prime
    q_n = q;
    q_prime_n = q_prime;
    % =====================================================================
    % ======= MAIN ROUTINE
    % =====================================================================
    % === PREDICTOR
    % f(t,x(j-1))
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = q_store(i,4);
            q_prime(i,1) = q_prime_store(i,4);  
            m=m+1;
        end
    end 
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,h);
    if er==1
        fprintf('error while calculating f1 for A-B-M\n');
        y_out=zeros(2*(7*nBody-nCon),1);
        return;        
    end
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f1 = [v_p_out;v_pp_out];    
    % f(t-h,x(j-2))
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = q_store(i,3);
            q_prime(i,1) = q_prime_store(i,3); 
            m=m+1;
        end
    end 
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t-h,h);
    if er==1
        fprintf('error while calculating f2 for A-B-M\n'); 
        y_out=zeros(2*(7*nBody-nCon),1);
        return;        
    end
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f2 = [v_p_out;v_pp_out];  
    % f(t-2*h,x(j-3))
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = q_store(i,2);
            q_prime(i,1) = q_prime_store(i,2);
            m=m+1;
        end
    end 
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t-2*h,h);
    if er==1
        fprintf('error while calculating f3 for A-B-M\n');     
        y_out=zeros(2*(7*nBody-nCon),1);
        return;        
    end
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f3 = [v_p_out;v_pp_out]; 
    % f(t-3*h, x(j-4))
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = q_store(i,1);
            q_prime(i,1) = q_prime_store(i,1); 
            m=m+1;
        end
    end 
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t-3*h,h);
    if er==1
        fprintf('error while calculating f4 for A-B-M\n'); 
        y_out=zeros(2*(7*nBody-nCon),1);
        return;        
    end
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f4 = [v_p_out;v_pp_out];
    % ============== PREDICTOR
    y_pred = y_n + (h/24)*(55*f1 - 59*f2 + 37*f3 - 9*f4);
    % one more evaluation
    % f(t+h, x_pred)
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q(i,1) = y_pred(m+1,1);
            q_prime(i,1) = y_pred(nind+m+1,1); 
            m=m+1;
        end
    end 
    [~,q_prime_out,q_prime_prime_out,~,~,~,~,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t+h,h);
    if er==1
        fprintf('error while calculating f5 for A-B-M\n'); 
        y_out=zeros(2*(7*nBody-nCon),1);
        return;        
    end
    m=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_p_out(m+1,1) = q_prime_out(i,1);
            v_pp_out(m+1,1) = q_prime_prime_out(i,1);
            m=m+1;
        end
    end
    f5 = [v_p_out;v_pp_out];
    % ============== CORRECTOR
    y_out = y_n + (h/24)*(9*f5 + 19*f1 - 5*f2 + f3);
    
%     xbar = x(j-1)+(h/24)*( 55*f(t,x(j-1)) - 59*f(t-h,x(j-2))+ 37*f(t-2*h,x(j-3)) - 9*f(t-3*h, x(j-4)) ); % predictor
%     x(j)=x(j-1) + (h/24)*( 9*f(t+h, xbar) + 19*f(t, x(j-1))-5*f(t-h, x(j-2)) +f(t-2*h, x(j-3)));% corrector
%     t = ta + (j-1)*h;


    
end