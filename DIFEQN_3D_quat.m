function [q_out,q_prime_out,q_prime_prime_out,FA,MA,MAT,lambda,er] = DIFEQN_3D_quat(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,dt)
% [q_out,q_prime_out,q_prime_prime_out,Q,lambda,er]
    % =====================================================================
    % ======= INITIALIZE VARIABLES
    % =====================================================================
    % === allocate output
%     q_out = zeros(7*nBody,1); <- not necessary
    q_prime_out = zeros(7*nBody,1);
    q_prime_prime_out = zeros(7*nBody,1);    
    % === calculate number of constrains for Phi_pi
    nConstraints = nCon - nBody;
    % =====================================================================
    % ======= solve Phi(u,v)=0 using N-R
    % =====================================================================
    [q_out,~,~,~,er] =  N_R_dynamic(q,q_prime,n, joint,nBody,nJoint, nCon, t,dt);
%     [q_out,er] = N_R_dynamic(q,n, joint,nBody,nJoint, nCon, t);
    if er==1 
       fprintf('error at time t=%f, \n',t);
       FA=0;
       MA=0;
       MAT=0;
       lambda=0;
       return; 
    end
    % =====================================================================
    % ======= solve Phi_u*u_prime = -Phi_v*v_prime - Phi_t for u_prime
    % =====================================================================
    % === number of independent variables
    nind = 7*nBody-nCon;
    % === allocate Phi_u
    Phi_u = zeros(nCon,nCon);
    % === allocate Phi_v
    Phi_v = zeros(nCon,nind);
    % === allocate u_prime and v_prime
    u_prime = zeros(nCon,1);
    v_prime = zeros(7*nBody-nCon,1);    
    % === calc Phi_q
    Phi_q = calc_Phi_q(joint,q_out,nBody,nJoint,nCon);
    % === extract Phi_u and Phi_v
    m=0;
    l=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            Phi_v(:,m+1) = Phi_q(:,i);
            m=m+1;
        else
            l=l+1;
            Phi_u(:,l) = Phi_q(:,i); 
        end
    end
    % === extract u_prime and v_prime
    m=0;
    l=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            v_prime(m+1,1) = q_prime(i,1);
            m=m+1;
        else
            l=l+1;
            u_prime(l,1) = q_prime(i,1); 
        end
    end
    % === calc Phi_t
    Phi_t = calc_Phi_t(joint,nJoint,nCon);
    % === calc RHS = -Phi_v*v_prime - Phi_t  to be solved
    RHS = -Phi_v*v_prime - Phi_t;
    % === construct matrix to solve
    MAT = Phi_u;
    for i=1:nCon
        MAT(i,nCon+1) = RHS(i,1);
    end
    % ================ SOLVE ================
    [A, loc_indx, er] = RowReducedEchelonForm(MAT,1);
    if er ==1
        return;
    end
    u_prime = get_solution(A, nCon,nCon+1,loc_indx); 
    % ======= return q_prime =======
    m=0;
    l=0;
    for i=1:7*nBody
        if (m<nind && i==n(m+1,1))
            q_prime_out(i,1) = v_prime(m+1,1);
            m=m+1;
        else
            l=l+1;
            q_prime_out(i,1) = u_prime(l,1);             
        end
    end
    clear MAT RHS
    % =====================================================================
    % ======= solve DYNAMICS
    % =====================================================================
    % allocate and calculate G matrices
    for i=1:nBody       
        for j=1:4
            p{i}(j,1)= q_out(j+(i-1)*7+3);
            p_prime{i}(j,1) = q_prime_out(j+(i-1)*7+3);
        end
    end
    for i=1:nBody
        G{i} = zeros(3,4); 
        G_p{i} = zeros(3,4);
        G{i} = calc_G(p{i});
        G_p{i} = calc_G(p_prime{i});
    end  
    % === calculate gamma vector (already has quaternions)
    gamma = calc_gamma(joint,q_out,q_prime_out,nBody,nJoint,nCon);
    % === calculate Phi_q
    Phi_q = calc_Phi_q(joint,q_out,nBody,nJoint,nCon);
    % === split Phi_q on Phi_r, Phi_p and Phi_pp
    for j=1:nBody
        m=(j-1)*3;
        l = (j-1)*7;
        for i=1:nConstraints
            Phi_r(i,m+1) = Phi_q(i,l+1); 
            Phi_r(i,m+2) = Phi_q(i,l+2); 
            Phi_r(i,m+3) = Phi_q(i,l+3);  
        end        
    end
    for j=1:nBody
        m=(j-1)*4;
        l = (j-1)*7+3;
        for i=1:nConstraints
            Phi_p(i,m+1) = Phi_q(i,l+1); 
            Phi_p(i,m+2) = Phi_q(i,l+2); 
            Phi_p(i,m+3) = Phi_q(i,l+3); 
            Phi_p(i,m+4) = Phi_q(i,l+4); 
        end        
    end
    Phi_pp = zeros(nBody,4*nBody);
    for i=1:nBody    
        m = (i-1)*4;
        Phi_pi = 2.0*p{i}';
        for j=1:4
            Phi_pp(i,m+j) = Phi_pi(j);
        end                           
    end
    % === assemble M matrix (mass)
    M = zeros(3*nBody,3*nBody);
    for i=1:nBody
        m=(i-1)*3;
        M(m+1,m+1) = body(i).m;
        M(m+2,m+2) = body(i).m;
        M(m+3,m+3) = body(i).m;
    end
    % === assemble 4*G'JG matrix (J is moment of inertia)
    GJG = zeros(4*nBody,4*nBody);
    for i=1:nBody
        l=(i-1)*4;        
        GJG(l+1:l+4,l+1:l+4) = 4*G{i}'*body(i).I*G{i};
    end
    % =====================================================================
    % ======= assemble LHS
    LHS = [M,zeros(3*nBody,4*nBody),Phi_r',zeros(3*nBody,nBody);
           zeros(4*nBody,3*nBody), GJG, Phi_p', Phi_pp';
           Phi_r,Phi_p,zeros(nConstraints,nConstraints),zeros(nConstraints,nBody);
           zeros(nBody,3*nBody),Phi_pp,zeros(nBody,nConstraints),zeros(nBody,nBody)];
    % =====================================================================
    % ========== assign Applied Forces vector ==========    
%     FA = calc_FA(q_out,q_prime_out,nBody,body);
    [FA,~,~,~,~,~] = calc_FA(q,q_prime,nBody,body);
    % =====================================================================
    % ========== assign Applied Moments vector ==========
%     MA = calc_MA(nBody,body);
    MA = calc_MA(q,q_prime,nBody,body);
    MOM = zeros(4*nBody,1);
    for i=1:nBody
        m=(i-1)*4;
        l=(i-1)*3;
        MOM(m+1:m+4,1) = 2*G{i}'*MA(l+1:l+3,1) + 8*G_p{i}'*body(i).I*G_p{i}*p{i};
    end
    % =====================================================================
    % === assemble RHS vector
    for i=1:3*nBody
        RHS(i,1) = FA(i,1);
    end
    m=0;
    for i=(3*nBody+1):7*nBody
        m=m+1;
        RHS(i,1) = MOM(m,1);
    end    
    m=0;
    for i=(7*nBody+1):(7*nBody+nCon)
        m=m+1;
        RHS(i,1) = gamma(m,1);    
    end
    % === assemble matrix to solve
    MAT = LHS;
    MAT(:,(7*nBody+nCon+1)) = RHS;
    % =====================================================================
    % ================== SOLVE ============================================
    [A, loc_indx, er] = RowReducedEchelonForm(MAT,1);
    SOL = get_solution(A, 7*nBody+nCon,7*nBody+nCon+1,loc_indx);
    % extract lambdas
    lambda(:,1) = SOL((7*nBody+1):(7*nBody+nConstraints),1);
    % extract accelerations
    q_prime_prime_sol = SOL(1:(7*nBody),1);
    % === extract q_prime_prime
    for i=1:nBody
        r_s = (i-1)*7;
        r_sol_s = (i-1)*3;
        q_prime_prime_out(r_s+1:r_s+3,1) = q_prime_prime_sol(r_sol_s+1:r_sol_s+3,1);
        q_prime_prime_out(r_s+4:r_s+7,1) = q_prime_prime_sol(3*nBody+r_sol_s+1:3*nBody+r_sol_s+4,1);
    end
    
   
    %%%%%%%%%%%%%%%%% TMP %%%%%%%%%%%%%%%%%%%%%%%%%%
%     FA=LHS;
%     MA=Phi_r;
%     MAT=RHS;
%     lambda=0;
    %%%%%%%%%%%%%%%%% TMP %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

%     % === assemble M matrix (mass)
%     M = zeros(3*nBody,3*nBody);
%     for i=1:nBody
%         m=(i-1)*3;
%         M(m+1,m+1) = body(i).m;
%         M(m+2,m+2) = body(i).m;
%         M(m+3,m+3) = body(i).m;
%     end
%     % === assemble J matrix (moments of inertia)
%     J = zeros(3*nBody,3*nBody);
%     for i=1:nBody
%         m = (i-1)*3;
%         J(m+1:m+3,m+1:m+3) = body(i).I;
%     end
%     % === extract Phi_r matrix (part of Jacobian with position derivatives)
%     Phi_r = zeros(nConstraints,3*nBody);
%     for i=1:nBody
%         m=(i-1)*3;
%         n = (i-1)*6;
%         Phi_r(:,m+1:m+3) = Phi_q_pi(:,n+1:n+3);
%     end
%     % === extract Phi_pi matrix (part of Jacobian with angle derivatives)
%     Phi_pi = zeros(nConstraints,3*nBody);
%     for i=1:nBody
%         m=(i-1)*3;
%         n = (i-1)*6;
%         Phi_pi(:,m+1:m+3) = Phi_q_pi(:,n+4:n+6);
%     end
%     % ============ assemble LHS ============
%     LHS = zeros(6*nBody+nConstraints,6*nBody+nConstraints);
%     for i=1:3*nBody
%         for j=1:3*nBody
%             LHS(i,j) = M(i,j);
%         end
%     end
%     m=0;
%     for i=(3*nBody+1):(6*nBody)
%         m=m+1;
%         n=0;
%         for j=(3*nBody+1):(6*nBody)
%             n=n+1;
%             LHS(i,j) = J(m,n);
%         end
%     end
%     m=0;
%     for i=(6*nBody+1):(6*nBody+nConstraints)
%         m=m+1;
%         n=0;
%         for j=1:(3*nBody)
%             n=n+1;
%             LHS(i,j) = Phi_r(m,n);
%         end
%     end
%     m=0;
%     for i=(6*nBody+1):(6*nBody+nConstraints)
%         m=m+1;
%         n=0;
%         for j=(3*nBody+1):(6*nBody)
%             n=n+1;
%             LHS(i,j) = Phi_pi(m,n);
%         end
%     end
%     m=0;
%     for i=1:(3*nBody)
%         m=m+1;
%         n=0;
%         for j=(6*nBody+1):(6*nBody+nConstraints)
%             n=n+1;
%             LHS(i,j) = Phi_r(n,m);
%         end
%     end
%     m=0;
%     for i=(3*nBody+1):(6*nBody)
%         m=m+1;
%         n=0;
%         for j=(6*nBody+1):(6*nBody+nConstraints)
%             n=n+1;
%             LHS(i,j) = Phi_pi(n,m);
%         end
%     end
%     % === assemble scew_omegas vector
%     sc_om = zeros(3*nBody,3*nBody);
%     for i=1:nBody
%         m=(i-1)*3;
%         tmp = scew_sym(omega{i});
%         sc_om(m+1:m+3,m+1:m+3) = tmp;
%     end
%     % === assemble omegas vector
%     om = zeros(3*nBody,1);
%     for i=1:nBody
%         m=(i-1)*3;
%         om(m+1:m+3,1) = omega{i};
%     end
%     % ========== asign Applied Forces vector ==========    
% %     for i=1:nBody
% %         m=(i-1)*3;
% %         FA(m+1:m+3,1) = [0;-body(i).m*9.81;0];        
% %     end
%     FA = calc_FA(q_out,q_prime_out,nBody,body);
%     % ========== asign Applied Moments vector ==========
% %     for i=1:3*nBody
% %         MA(i,1) = 0.0;
% %     end
%     MA = calc_MA(nBody);
%     MOM = MA - sc_om*J*om;
% %     MOM
% %     sc_om
% %     om
% %     J
%     % === assemble RHS vector
%     for i=1:3*nBody
%         RHS(i,1) = FA(i,1);
%     end
%     m=0;
%     for i=(3*nBody+1):6*nBody
%         m=m+1;
%         RHS(i,1) = MOM(m,1);
%     end    
%     m=0;
%     for i=(6*nBody+1):(6*nBody+nConstraints)
%         m=m+1;
%         RHS(i,1) = gamma(m,1);    
%     end
%     % === assemble matrix to solve
%     MAT = LHS;
%     MAT(:,(6*nBody+nConstraints+1)) = RHS;
%     % === SOLVE
%     [A, loc_indx, er] = RowReducedEchelonForm(MAT,1);
%     SOL = get_solution(A, 6*nBody+nConstraints,6*nBody+nConstraints+1,loc_indx);
%     if er ==1
%         return;
%     end
%     % extract lambdas
%     lambda(:,1) = SOL((6*nBody+1):(6*nBody+nConstraints),1);
%     % ========== extract q_prime_prime_out ==========
%     q_prime_prime_PI = SOL(1:(6*nBody),1);
%     % === extract omega_prime
%     for i=1:nBody
%         m=3*nBody + (i-1)*3;
%         omega_prime{i} = q_prime_prime_PI(m+1:m+3,1);
%     end
%     
%     for i=1:nBody
%         r_sart = (i-1)*7;
%         r_pi_start = (i-1)*3;
%         q_prime_prime_out(r_sart+1:r_sart+3,1) = q_prime_prime_PI(r_pi_start+1:r_pi_start+3,1);
%         q_prime_prime_out(r_sart+4:r_sart+7,1) = 0.5*G{i}'*omega_prime{i} - 0.25*omega{i}'*omega{i}*p{i};
%     end
%     q_prime_prime_out = q_prime_prime_PI;
%     omega_prime{1}
%     omega_prime{2}
%     omega_prime{3}
%     omega_prime{4}
    
    
    
    
end