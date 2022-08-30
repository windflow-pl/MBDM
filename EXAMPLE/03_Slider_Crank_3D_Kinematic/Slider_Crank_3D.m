clc
clear all

%% INSTALL MBDM FUNCTIONS
addpath('../../.');

%% input
nBody = 4;
nJoint = 7; % include driving

% =================================
% ======== allocate bodies ========
% =================================
for(i=1:nBody)
    body(i) = struct('name',{''},'r',{[0,0,0]},'p',{[0,0,0,0]},'v',{[0,0,0]},'omega',{[0,0,0]},'m',{0},'I',{zeros(3)},'A',{zeros(3)});
end
for(i=1:nJoint)
    joint(i) = struct('type',{''},'bodies',{[0,0]},'p1',{[0,0,0]},'q1',{[0,0,0]},'r1',{[0,0,0]},'p2',{[0,0,0]},'q2',{[0,0,0]},'r2',{[0,0,0]},'dist',{0});
end
% ===============================
% ======== define bodies ========
% ===============================
% === inital orientation
angle = [pi, 0, 0];%euler angles 1-2-3 notation
A = calc_A_angle(angle);
e0 = sqrt((trace(A) +1)/4);
if (e0 ~= 0)
    e1 = (A(3,2)-A(2,3))/(4*e0);
    e2 = (A(1,3)-A(3,1))/(4*e0);
    e3 = (A(2,1)-A(1,2))/(4*e0);
end
if (e0 == 0)
    e1 = sqrt((1+2*A(1,1)-trace(A))/4);
    e2 = sqrt((1+2*A(2,2)-trace(A))/4);
    e3 = sqrt((1+2*A(3,3)-trace(A))/4);
end
% === create body
body(1).name = 'crank';
body(1).r = [0.0,0.1,0.12];
body(1).p = [e0,e1,e2,e3];
% body(1).p = [0.7,0.71,0.0,0.0];
body(1).v = [0,0,0];
body(1).omega = [0,0,2*pi];
body(1).m = 0.12;
body(1).I = [0.0001,0,0; 0,0.00001,0; 0,0,0.0001];
body(1).A = calc_A(body(1).p');
% === inital orientation
% angle = [0.0, 0.0, -0.4];
% A = calc_A_angle(angle);
% e0 = sqrt((trace(A) +1)/4);
% if (e0 ~= 0)
%     e1 = (A(3,2)-A(2,3))/(4*e0);
%     e2 = (A(1,3)-A(3,1))/(4*e0);
%     e3 = (A(2,1)-A(1,2))/(4*e0);
% end
% if (e0 == 0)
%     e1 = sqrt((1+2*A(1,1)-trace(A))/4);
%     e2 = sqrt((1+2*A(2,2)-trace(A))/4);
%     e3 = sqrt((1+2*A(3,3)-trace(A))/4);
% end
% === create body
body(2).name = 'rod';
body(2).r = [0.1,0.05,0.1];
% body(2).p = [e0,e1,e2,e3];
body(2).p = [0.8,-0.21,0.4,-0.1];
body(2).v = [0,0,0];
body(2).omega = [0,0,0];
body(2).m = 0.5;
body(2).I = [0.004,0,0; 0,0.0004,0; 0,0,0.004];
body(2).A = calc_A(body(2).p');
% === inital orientation
angle = [0.0, 0.0, 0.0];
A = calc_A_angle(angle);
e0 = sqrt((trace(A) +1)/4);
if (e0 ~= 0)
    e1 = (A(3,2)-A(2,3))/(4*e0);
    e2 = (A(1,3)-A(3,1))/(4*e0);
    e3 = (A(2,1)-A(1,2))/(4*e0);
end
if (e0 == 0)
    e1 = sqrt((1+2*A(1,1)-trace(A))/4);
    e2 = sqrt((1+2*A(2,2)-trace(A))/4);
    e3 = sqrt((1+2*A(3,3)-trace(A))/4);
end
% === create body
body(3).name = 'slider';
body(3).r = [0.2,0.0,0.0];
body(3).p = [e0,e1,e2,e3];
body(3).v = [0,0,0];
body(3).omega = [0,0,0];
body(3).m = 2.0;
body(3).I = [0.0001,0,0; 0,0.0001,0; 0,0,0.0001];
body(3).A = calc_A(body(3).p');
% === inital orientation
angle = [0.0, 0.0, 0.0];
A = calc_A_angle(angle);
e0 = sqrt((trace(A) +1)/4);
if (e0 ~= 0)
    e1 = (A(3,2)-A(2,3))/(4*e0);
    e2 = (A(1,3)-A(3,1))/(4*e0);
    e3 = (A(2,1)-A(1,2))/(4*e0);
end
if (e0 == 0)
    e1 = sqrt((1+2*A(1,1)-trace(A))/4);
    e2 = sqrt((1+2*A(2,2)-trace(A))/4);
    e3 = sqrt((1+2*A(3,3)-trace(A))/4);
end
% === create body
body(4).name = 'ground';
body(4).r = [0.0,0.0,0.0];
body(4).p = [e0,e1,e2,e3];
body(4).v = [0,0,0];
body(4).omega = [0,0,0];
body(4).m = 1;
body(4).I = [1,0,0; 0,1,0; 0,0,1];
body(4).A = calc_A(body(4).p');
% ====== clear
clear angle A e0 e1 e2 e3
% ===============================
% ======== define joints ========
% ===============================
% % revolute joint
joint(1).type = 'R';
joint(1).bodies = [1,4];%IDs of body 1 and body 2
joint(1).p1 = [0,0,0];% joint location in LCS of body 1
joint(1).q1 = [1,0,0];% point on the axis of rotation in LCS 1
joint(1).r1 = [0,1,0];% point located on the normal to the axis of rotation in LCS 1
joint(1).p2 = [0,0.1,0.12];% joint location in LCS of body 2
joint(1).q2 = [1,0.1,0.12];% point on the axis of rotation in LCS 2
% % spherical joint
joint(2).type = 'S';
joint(2).bodies = [1,2];%IDs of body 1 and body 2
joint(2).p1 = [0,0.08,0];% joint location in LCS of body 1
joint(2).p2 = [-0.15,0,0];% joint location in LCS of body 2
% % revolute-cylindrical joint
joint(3).type = 'RC';
joint(3).bodies = [2,3];% body 1 has a revolute part and body 2 translational part
joint(3).p1 = [0.15,0,0];% joint location in LCS of body 1
joint(3).q1 = [0.15,1,0];% point on the axis of rotation in LCS 1
joint(3).r1 = [1.15,0,0];% point located on the normal to the axis of rotation in LCS 1
joint(3).p2 = [0.2,0,0];% any point on the axis of translation in LCS of body 2
joint(3).q2 = [1.2,0,0];% offset point on the axis of translation in LCS 2
joint(3).r2 = [0.2,1,0];% point located on the normal to the axis of translation in LCS 1
% % translational joint
joint(4).type = 'T';
joint(4).bodies = [3,4];%IDs of body 1 and body 2
joint(4).p1 = [0,0,0];% any point in the axis of translation in LCS of body 1
joint(4).q1 = [1,0,0];% offset point in the axis of translation in LCS 1
joint(4).r1 = [0,1,0];% offset point normal to the axis of translation in LCS 1
joint(4).p2 = [0.2,0,0];% any point in the axis of translation in LCS of body 2
joint(4).q2 = [1.2,0,0];% offset point in the axis of translation in LCS 2
joint(4).r2 = [0.2,0,1];% offset point normal to the axis of translation and normal to r1 in LCS 2
% % disntance constraint
joint(5).type = 'DI';
joint(5).bodies = [2,3];%IDs of body 1 and body 2
joint(5).p1 = [0.15,0,0];% point in LCS of body 1
joint(5).p2 = [1.0,0,0];% point in LCS of body 2
joint(5).dist = 1.0; % desired distance netween those two points
% % ground
joint(6).type = 'G';
joint(6).bodies = [4,0];%ID of the body to be fixed to the ground
joint(6).p1 = [0,0,0];% desired location of LCS in GCS
joint(6).q1 = [0,0,0];% desired orientation in terms of Euler parameters [p2,p3,p4]
% % Driving
joint(7).type = 'RDRV';% rotational driving
joint(7).bodies = [4,1];%ID 2 is ID of the body to be driven relative to body 1
joint(7).p2 = [0,0,1];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
% angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
joint(7).p1 = [0,1,0];% versor in LCS 1 that lies in plane of rotation (see above)
joint(7).q1 = [1,0,0];% versor in LCS 1 indicating the axis of rotation
joint(7).r1 = [pi,0,0];% inital angle in the first entry (from p1 to p2) 0<angle<2*pi
joint(7).dist = -2*pi;% rotational velocity (positive counter-clockwise)

%% ======== calculate DOFs  ========
% ==================================
DOFs = nBody*6;
nConstraints = 0;
for(i=1:(nJoint))
    if(strcmp(joint(i).type,'R')) nConstraints =nConstraints+5; end
    if(strcmp(joint(i).type,'S')) nConstraints =nConstraints+3; end
    if(strcmp(joint(i).type,'RC')) nConstraints =nConstraints+3; end
    if(strcmp(joint(i).type,'T')) nConstraints =nConstraints+5; end
    if(strcmp(joint(i).type,'DI')) nConstraints =nConstraints+1; end
    if(strcmp(joint(i).type,'G')) nConstraints =nConstraints+6; end
    if(strcmp(joint(i).type,'RDRV')) nConstraints =nConstraints+1; end
end
disp(sprintf('Number of bodies = %d, Joints + Driving = %d',nBody,nJoint));
disp(sprintf('Number of constraint equations = %d, DOFs = %d',nConstraints,DOFs-nConstraints));
% == number of constraints with quaternions
nCon = nConstraints+nBody;

%% =========    KINEMATICS     ==========
% =======================================
t = 0.01;
% === allocate q
q = zeros(nBody*7,1);
% === create q
% q = [r1,p1,r2,p2,r3,p3...]
for i=1:nBody
    for j=1:3
        q(j+(i-1)*7,1) = body(i).r(j);
    end
    for j=1:4
        q(j+(i-1)*7+3,1) = body(i).p(j);
    end
end



Phi = calc_Phi(joint,q,t,nBody,nJoint,nCon);
Phi_q = calc_Phi_q(joint,q,nBody,nJoint,nCon);



% =======================================
% ======= define time step
% =======================================
dt = 0.01;
t_final = 1;
% =======================================
% ======= define tracers
% =======================================
t11 = [0.0, 0.0 ,0.0];
t12 = [0.0, 0.08, 0.0];
t21 = [-0.15, 0, 0];
t22 = [0.15, 0, 0];
t31 = [0.02, -0.01, -0.01];
t32 = [0.02, -0.01, 0.01];
t33 = [0.02, 0.01, 0.01];
t34 = [0.02, 0.01, -0.01];
t35 = [-0.02, -0.01, -0.01];
t36 = [-0.02, -0.01, 0.01];
t37 = [-0.02, 0.01, 0.01];
t38 = [-0.02, 0.01, -0.01];

t41 = [0.0,0.1,0.12]+[0.01, -0.01, -0.01];
t42 = [0.0,0.1,0.12]+[0.01, -0.01, 0.01];
t43 = [0.0,0.1,0.12]+[0.01, 0.01, 0.01];
t44 = [0.0,0.1,0.12]+[0.01, 0.01, -0.01];
t45 = [0.0,0.1,0.12]+[-0.01, -0.01, -0.01];
t46 = [0.0,0.1,0.12]+[-0.01, -0.01, 0.01];
t47 = [0.0,0.1,0.12]+[-0.01, 0.01, 0.01];
t48 = [0.0,0.1,0.12]+[-0.01, 0.01, -0.01];

% =======================================
% ======= MAIN LOOP
% =======================================
for k=1:(t_final/dt+1)
    t = (k-1)*dt;
    % ======== clear and allocate q
    q = zeros(nBody*7,1);
    % ======== get q
    for i=1:nBody
        for j=1:3
            q(j+(i-1)*7,1) = body(i).r(j);
        end
        for j=1:4
            q(j+(i-1)*7+3,1) = body(i).p(j);
        end
    end
    % ======== calculate q
    [q_out,er] = N_R_kinematic(q, joint,nBody,nJoint, nCon, t);
    if er==1
        return;
    end
    % ======== calculate q_prime
    q_prime = calc_q_p(q_out,joint,nBody,nJoint,nCon);
    % ======== calculate q_prime_prime
    q_prime_prime = calc_q_pp(q_out,q_prime,joint,nBody,nJoint,nCon);
    % ======== update Bodies
    for i=1:nBody
        for j=1:3
            body(i).r(j) = q_out(j+(i-1)*7,1);
            body(i).v(j) = q_prime(j+(i-1)*7,1);
        end
        for j=1:4
            body(i).p(j) = q_out(j+(i-1)*7+3,1);
            p_prime{i}(j,1) = q_prime(j+(i-1)*7+3,1);
        end
        body(i).A = calc_A(body(i).p');
        body(i).omega = (2*calc_G(body(i).p')*p_prime{i})';
    end
    % === create output of slider variables
    x(k,1) = t;
    x(k,2) = q_out(15,1);
    x(k,3) = q_prime(15,1);
    x(k,4) = q_prime_prime(15,1);
    % =====================================================================
    % ======== END OF KINEMATICS ==========================================
    % ======== START INVERT DYNAMICS ======================================
    % =====================================================================
    % === calculate gamma vector
    gamma = calc_gamma(joint,q_out,q_prime,nBody,nJoint,nCon);
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
        tmp = scew_sym(body(i).omega');
        sc_om(m+1:m+3,m+1:m+3) = tmp;
    end
    % === assemble omegas vector
    om = zeros(3*nBody,1);
    for i=1:nBody
        m=(i-1)*3;
        om(m+1:m+3,1) = body(i).omega';
    end
    % ========== asign Applied Forces vector ==========
    m=0;
    m2 = 0;
    for i=1:3*nBody
        m=m+1;
        if m==3
            m=0;
            m2 = m2+1;
            FA(i,1) = -body(m2).m*9.81;
        else
            FA(i,1) = 0.0;
        end
    end
    % ========== asign Applied Moments vector ==========
    for i=1:3*nBody
        MA(i,1) = 0.0;
    end
    MOM = MA - sc_om*J*om;
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
    % extract lambdas
    lambda(:,1) = SOL((6*nBody+1):(6*nBody+nConstraints),1);
    % =====================================================================
    % ======== END OF INVERT DYNAMICS =====================================
    % =====================================================================
    % ======== EXTRACT FORCES AND MOMENT ON CRANK JOINT
##    C_1_DRV = eye(3);
##    hi = (joint(1).q1 - joint(1).p1)';
##    fi = (joint(1).r1 - joint(1).p1)';
##    gi = -(scew_sym(fi)*hi);
##    C_1_REV = [fi,gi,hi];
##    lambda_rev = lambda(1:5,1);
##    lambda_drv = lambda(24,1);
##    Phi_ri_rev = Phi_r(1:5,1:3);
##    Phi_ri_drv = Phi_r(24,1:3);
##    Phi_pi_drv = Phi_pi(24,1:3);
##    % === calculate forces
##%     F = body(1).A*C_1_REV*(-C_1_REV'*body(1).A'*Phi_ri_rev'*lambda_rev);
##%     T = body(1).A*C_1_REV*(-C_1_DRV'*(Phi_pi_drv' - scew_sym([0;0;0])*body(1).A'*Phi_ri_drv')*lambda_drv);
##    F = (-C_1_REV'*body(1).A'*Phi_ri_rev'*lambda_rev);
##    T = (-C_1_DRV'*(Phi_pi_drv' - scew_sym([0;0;0])*body(1).A'*Phi_ri_drv')*lambda_drv);
##    force(k,1) = t;
##    force(k,2) = F(1,1);
##    force(k,3) = F(2,1);
##    force(k,4) = F(3,1);
##    force(k,5) = T(1,1);
##    force(k,6) = T(2,1);
##    force(k,7) = T(3,1);


    % ===============================
    % ======== trace and plot tracers
    % ===============================
    Tr11 = body(1).r' + body(1).A*t11';
    Tr12 = body(1).r' + body(1).A*t12';
    Tr21 = body(2).r' + body(2).A*t21';
    Tr22 = body(2).r' + body(2).A*t22';
    Tr31 = body(3).r' + body(3).A*t31';
    Tr32 = body(3).r' + body(3).A*t32';
    Tr33 = body(3).r' + body(3).A*t33';
    Tr34 = body(3).r' + body(3).A*t34';
    Tr35 = body(3).r' + body(3).A*t35';
    Tr36 = body(3).r' + body(3).A*t36';
    Tr37 = body(3).r' + body(3).A*t37';
    Tr38 = body(3).r' + body(3).A*t38';

    Tr41 = body(4).r' + body(4).A*t41';
    Tr42 = body(4).r' + body(4).A*t42';
    Tr43 = body(4).r' + body(4).A*t43';
    Tr44 = body(4).r' + body(4).A*t44';
    Tr45 = body(4).r' + body(4).A*t45';
    Tr46 = body(4).r' + body(4).A*t46';
    Tr47 = body(4).r' + body(4).A*t47';
    Tr48 = body(4).r' + body(4).A*t48';

    crank(1,:) = Tr11;
    crank(2,:) = Tr12;
    arm(1,:) = Tr21;
    arm(2,:) = Tr22;
    slider(1,:) = Tr31;
    slider(2,:) = Tr32;
    slider(3,:) = Tr33;
    slider(4,:) = Tr34;
    slider(5,:) = Tr35;
    slider(6,:) = Tr36;
    slider(7,:) = Tr37;
    slider(8,:) = Tr38;

    ground(1,:) = Tr41;
    ground(2,:) = Tr42;
    ground(3,:) = Tr43;
    ground(4,:) = Tr44;
    ground(5,:) = Tr45;
    ground(6,:) = Tr46;
    ground(7,:) = Tr47;
    ground(8,:) = Tr48;

     hold on
     clf

     grid on
##     plot3(crank(:,1),crank(:,2),crank(:,3),'o-');
     plot3(crank(:,1),crank(:,2),crank(:,3),'o-',arm(:,1),arm(:,2),arm(:,3),'s-',slider(:,1),slider(:,2),slider(:,3),'+-',ground(:,1),ground(:,2),ground(:,3),'+-' );
     axis([-0.05 0.4 -0.2 0.2 -0.05 0.25])
      daspect([1 1 1])
     xlabel('x')
     ylabel('y')
     zlabel('z')
     title(['t=', num2str(t,'%1.5f\n')])
##     view(115,20)

     hold off
     drawnow
##     pause(0.001);



end







%% =========== PLOTS

% %% plot against the book
src1 = 'results/data/x.dat';
src2 = 'results/data/x_p.dat';
src3 = 'results/data/x_pp.dat';
delimiterIn = ' ';
x_book = load(src1);
x_p_book = load(src2);
x_pp_book = load(src3);

figure
hold on
plot(x(:,1),x(:,2),'k',x_book(:,1),x_book(:,2),'r');
hold off
figure
hold on
plot(x(:,1),x(:,3),'k',x_p_book(:,1),x_p_book(:,2),'r');
hold off
figure
hold on
plot(x(:,1),x(:,4),'k',x_pp_book(:,1),x_pp_book(:,2),'r');
hold off
%% save
##src = 'results/data/x_xp_xpp_slider.dat';
##save(src, 'x', '-ascii');






