function [q_out,q_prime_out,q_prime_prime_out,FA,MA,MAT,lambda,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,dt)
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
    % === allocate Phi_u
    Phi_u = zeros(nCon,nCon);
    % === allocate Phi_v
    Phi_v = zeros(nCon,7*nBody-nCon);
    % === allocate u_prime and v_prime
    u_prime = zeros(nCon,1);
    v_prime = zeros(7*nBody-nCon,1);
    % === number of independent variables
    nind = 7*nBody-nCon;
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
    % === calculate omegas
    % alocate p and p_prime
    for i=1:nBody
        p{i} = zeros(4,1);
        p_prime{i} = zeros(4,1);
    end
    % strip q_out and q_prime_out
    for i=1:nBody
        for j=1:4
            p{i}(j,1)= q_out(j+(i-1)*7+3);
            p_prime{i}(j,1)= q_prime_out(j+(i-1)*7+3);
        end
    end
    % allocate G matrices
    for i=1:nBody
        G{i} = zeros(3,3);
        G{i} = calc_G(p{i});
    end
    for i=1:nBody
        omega{i} = 2.0*G{i}*p_prime{i};
    end
    % === calculate gamma vector
    gamma = calc_gamma(joint,q_out,q_prime_out,nBody,nJoint,nCon);
    % gamma has quaternion constrains, remove
    gamma2 = gamma;
    clear gamma
    for i=1:nConstraints
        gamma(i,1) = gamma2(i,1);
    end
    clear gamma2
    % === calculate Phi_q_pi - Jacobian in angle variables, not quaternions
    Phi_q_pi = calc_Phi_q_pi(joint,q_out,nBody,nJoint,nConstraints);
    % === assemble M matrix (mass)
    M = zeros(3*nBody,3*nBody);
    for i=1:nBody
        m=(i-1)*3;
        M(m+1,m+1) = body(i).m;
        M(m+2,m+2) = body(i).m;
        M(m+3,m+3) = body(i).m;
    end
    % === assemble J matrix (moments of inertia)
    J = zeros(3*nBody,3*nBody);
    for i=1:nBody
        m = (i-1)*3;
        J(m+1:m+3,m+1:m+3) = body(i).I;
    end
    % === extract Phi_r matrix (part of Jacobian with position derivatives)
    Phi_r = zeros(nConstraints,3*nBody);
    for i=1:nBody
        m=(i-1)*3;
        n = (i-1)*6;
        Phi_r(:,m+1:m+3) = Phi_q_pi(:,n+1:n+3);
    end
    % === extract Phi_pi matrix (part of Jacobian with angle derivatives)
    Phi_pi = zeros(nConstraints,3*nBody);
    for i=1:nBody
        m=(i-1)*3;
        n = (i-1)*6;
        Phi_pi(:,m+1:m+3) = Phi_q_pi(:,n+4:n+6);
    end
    % ============ assemble LHS ============
    LHS = zeros(6*nBody+nConstraints,6*nBody+nConstraints);
    for i=1:3*nBody
        for j=1:3*nBody
            LHS(i,j) = M(i,j);
        end
    end
    m=0;
    for i=(3*nBody+1):(6*nBody)
        m=m+1;
        n=0;
        for j=(3*nBody+1):(6*nBody)
            n=n+1;
            LHS(i,j) = J(m,n);
        end
    end
    m=0;
    for i=(6*nBody+1):(6*nBody+nConstraints)
        m=m+1;
        n=0;
        for j=1:(3*nBody)
            n=n+1;
            LHS(i,j) = Phi_r(m,n);
        end
    end
    m=0;
    for i=(6*nBody+1):(6*nBody+nConstraints)
        m=m+1;
        n=0;
        for j=(3*nBody+1):(6*nBody)
            n=n+1;
            LHS(i,j) = Phi_pi(m,n);
        end
    end
    m=0;
    for i=1:(3*nBody)
        m=m+1;
        n=0;
        for j=(6*nBody+1):(6*nBody+nConstraints)
            n=n+1;
            LHS(i,j) = Phi_r(n,m);
        end
    end
    m=0;
    for i=(3*nBody+1):(6*nBody)
        m=m+1;
        n=0;
        for j=(6*nBody+1):(6*nBody+nConstraints)
            n=n+1;
            LHS(i,j) = Phi_pi(n,m);
        end
    end
    % === assemble scew_omegas vector
    sc_om = zeros(3*nBody,3*nBody);
    for i=1:nBody
        m=(i-1)*3;
        tmp = scew_sym(omega{i});
        sc_om(m+1:m+3,m+1:m+3) = tmp;
    end
    % === assemble omegas vector
    om = zeros(3*nBody,1);
    for i=1:nBody
        m=(i-1)*3;
        om(m+1:m+3,1) = omega{i};
    end
    % ========== asign Applied Forces vector ==========
%     for i=1:nBody
%         m=(i-1)*3;
%         FA(m+1:m+3,1) = [0;-body(i).m*9.81;0];
%     end
    FA = calc_FA(q_out,q_prime_out,nBody,body);
    % ========== asign Applied Moments vector ==========
%     for i=1:3*nBody
%         MA(i,1) = 0.0;
%     end
    MA = calc_MA(nBody);
    MOM = MA - sc_om*J*om;
%     MOM
%     sc_om
%     om
%     J
    % === assemble RHS vector
    for i=1:3*nBody
        RHS(i,1) = FA(i,1);
    end
    m=0;
    for i=(3*nBody+1):6*nBody
        m=m+1;
        RHS(i,1) = MOM(m,1);
    end
    m=0;
    for i=(6*nBody+1):(6*nBody+nConstraints)
        m=m+1;
        RHS(i,1) = gamma(m,1);
    end
    % === assemble matrix to solve
    MAT = LHS;
    MAT(:,(6*nBody+nConstraints+1)) = RHS;
    % === SOLVE
    [A, loc_indx, er] = RowReducedEchelonForm(MAT,1);
    SOL = get_solution(A, 6*nBody+nConstraints,6*nBody+nConstraints+1,loc_indx);
    if er ==1
        return;
    end
    % extract lambdas
    lambda(:,1) = SOL((6*nBody+1):(6*nBody+nConstraints),1);
    % ========== extract q_prime_prime_out ==========
    q_prime_prime_PI = SOL(1:(6*nBody),1);
    % === extract omega_prime
    for i=1:nBody
        m=3*nBody + (i-1)*3;
        omega_prime{i} = q_prime_prime_PI(m+1:m+3,1);
    end

    for i=1:nBody
        r_sart = (i-1)*7;
        r_pi_start = (i-1)*3;
        q_prime_prime_out(r_sart+1:r_sart+3,1) = q_prime_prime_PI(r_pi_start+1:r_pi_start+3,1);
        q_prime_prime_out(r_sart+4:r_sart+7,1) = 0.5*G{i}'*omega_prime{i} - 0.25*omega{i}'*omega{i}*p{i};
    end
%     q_prime_prime_out = q_prime_prime_PI;
%     omega_prime{1}
%     omega_prime{2}
%     omega_prime{3}
%     omega_prime{4}




end
