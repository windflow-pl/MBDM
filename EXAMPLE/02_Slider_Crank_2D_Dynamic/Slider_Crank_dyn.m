clc
clear all

%% INSTALL MBDM FUNCTIONS
addpath('../../.')

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
angle = [0, 0, pi];%euler angles in 1-2-3 notation (pitch-roll-yaw) (x-y-z)
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
body(1).r = [0.0,0.0,0.0];
body(1).p = [e0,e1,e2,e3];
% body(1).p = [0.7,0.71,0.0,0.0];
body(1).v = [0,0,0];
body(1).omega = [0,0,30];
body(1).m = 200;
body(1).I = [450,0,0; 0,450,0; 0,0,450];
body(1).A = calc_A(body(1).p');
% === inital orientation
% angle = [0.0, 0.0, -0.4158];
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
body(2).r = [-0.25,0.0,0.0];
% body(2).p = [e0,e1,e2,e3];
body(2).p = [1,0,0,0];
body(2).v = [0, -30, 0];
body(2).omega = [0,	0, 17.1429];
body(2).m = 35;
body(2).I = [35,0,0; 0,35,0; 0,0,35];
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
body(3).r = [1.5,0.0,0.0];
body(3).p = [e0,e1,e2,e3];
body(3).v = [0,	0,	0];
body(3).omega = [0,0,0];
body(3).m = 25;
body(3).I = [0.02,0,0; 0,0.02,0; 0,0,0.02];
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
joint(1).q1 = [0,0,1];% point on the axis of rotation in LCS 1
joint(1).r1 = [1,0,0];% point located on the normal to the axis of rotation in LCS 1
joint(1).p2 = [0,0,0];% joint location in LCS of body 2
joint(1).q2 = [0,0,1];% point on the axis of rotation in LCS 2
% % revolute joint
joint(2).type = 'R';
joint(2).bodies = [1,2];%IDs of body 1 and body 2
joint(2).p1 = [2,0,0];% joint location in LCS of body 1
joint(2).q1 = [2,0,1];% point on the axis of rotation in LCS 1
joint(2).r1 = [3,0,0];% point located on the normal to the axis of rotation in LCS 1
joint(2).p2 = [-1.75,0,0];% joint location in LCS of body 2
joint(2).q2 = [-1.75,0,1];% point on the axis of rotation in LCS 2
% % cylindrical joint
joint(3).type = 'C';
joint(3).bodies = [2,3];%IDs of body 1 and body 2
joint(3).p1 = [1.75,0,-1];% any point on the axis of rotation in LCS of body 1
joint(3).q1 = [1.75,0,0];% offset point on the axis of rotation in LCS 1
joint(3).r1 = [2.75,0,-1];% offset point normal to the axis of rotation in LCS 1
joint(3).p2 = [0,0,0];% any point on the axis of rotation in LCS of body 2
joint(3).q2 = [0,0,1];% offset point on the axis of rotation in LCS 2
joint(3).r2 = [1,0,0];% offset point normal to the axis of rotation and normal to r1 in LCS 2
% % absolute constraint on Y
joint(4).type = 'ABY';
joint(4).bodies = [3,0];%ID of the body
joint(4).dist = 0;% value of y to fix LCS
% % % absolute constraint on Z
joint(5).type = 'ABZ';
joint(5).bodies = [3,0];%ID of the body
joint(5).dist = 0;% value of z to fix LCS
% % ground
joint(6).type = 'G';
joint(6).bodies = [4,0];%ID of the body to be fixed to the ground
joint(6).p1 = [0,0,0];% desired location of LCS in GCS
joint(6).q1 = [0,0,0];% desired orientation in terms of Euler parameters [p2,p3,p4]
% % Driver as angle fix
joint(7).type = 'RDRV';% rotational driving
joint(7).bodies = [4,3];%ID 2 is ID of the body to be driven relative to body 1
joint(7).p2 = [1,0,0];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
% angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
joint(7).p1 = [1,0,0];% versor in LCS 1 that lies in plane of rotation (see above)
joint(7).q1 = [0,0,1];% versor in LCS 1 indicating the axis of rotation
joint(7).r1 = [0,0,0];% inital angle in the first entry (from p1 to p2) 0<angle<2*pi
joint(7).dist = 0;% rotational velocity (positive counter-clockwise)

%% ======== calculate DOFs  ========
% ==================================
DOFs = nBody*6;
nConstraints = 0;
for(i=1:(nJoint))
    if(strcmp(joint(i).type,'R')) nConstraints =nConstraints+5; end
    if(strcmp(joint(i).type,'S')) nConstraints =nConstraints+3; end
    if(strcmp(joint(i).type,'RC')) nConstraints =nConstraints+3; end
    if(strcmp(joint(i).type,'C')) nConstraints =nConstraints+4; end
    if(strcmp(joint(i).type,'T')) nConstraints =nConstraints+5; end
    if(strcmp(joint(i).type,'DI')) nConstraints =nConstraints+1; end
    if(strcmp(joint(i).type,'ABX')) nConstraints =nConstraints+1; end
    if(strcmp(joint(i).type,'ABY')) nConstraints =nConstraints+1; end
    if(strcmp(joint(i).type,'ABZ')) nConstraints =nConstraints+1; end
    if(strcmp(joint(i).type,'G')) nConstraints =nConstraints+6; end
    if(strcmp(joint(i).type,'RDRV')) nConstraints =nConstraints+1; end
end
disp(sprintf('Number of bodies = %d, Joints + Driving = %d',nBody,nJoint));
disp(sprintf('Number of constraint equations = %d, DOFs = %d',nConstraints,DOFs-nConstraints));
% == number of constraints with quaternions
nCon = nConstraints+nBody;

%% =========    KINEMATICS     ==========
% =======================================
t = 0.0;
% dt = 0.001;
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
% === create q_prime
q_prime = zeros(nBody*7,1);
for i=1:nBody
    for j=1:3
        q_prime(j+(i-1)*7,1) = body(i).v(j);
    end
    p_prime = 0.5*(calc_G(body(i).p'))'*body(i).omega';
    for j=1:4
        q_prime(j+(i-1)*7+3,1) = p_prime(j,1);
    end
end
clear p_prime
Phi = calc_Phi(joint,q,t,nBody,nJoint,nCon);
Phi_q = calc_Phi_q(joint,q,nBody,nJoint,nCon);
[A, loc_indx, er] = RowReducedEchelonForm(Phi_q,0);
% find independent variables
m=0;
for i=(nCon+1):7*nBody
    m=m+1;
    n(m,1) =  loc_indx(1,i);
end
fprintf('Independent variable assigned to n=%d\n',n(1,1));
% n = 7;
% number of independent variables
nind = 7*nBody-nCon;
% store variables for ABM scheme
q_store(:,1) = q;
q_prime_store(:,1) = q_prime;




% =======================================
% ======= define time step
% =======================================
dt = 0.005;
t_final =1;
% =======================================
% ======= define tracers
% =======================================
t11 = [0 , 0, 0];
t12 = [2 , 0, 0];
t21 = [-1.75, 0, 0];
t22 = [1.75, 0, 0];
t31 = [-0.2, -0.1, 0];
t32 = [-0.2, 0.1, 0];
t33 = [0.2, 0.1, 0];
t34 = [0.2, -0.1, 0];
t41 = [-0.1, -0.1, 0];
t42 = [-0.1, 0.1, 0];
t43 = [0.1, 0.1, 0];
t44 = [0.1, -0.1, 0];
% =======================================
% ======= MAIN LOOP (NUMERICAL INTEGRATION ROUTINE)
% =======================================
tmp = 0;
% n=0;
for k=1:(t_final/dt+1)
    t = (k-1)*dt;
    disp(sprintf('t = %f',t));
    % find independent variables
    Phi = calc_Phi(joint,q,t,nBody,nJoint,nCon);
    Phi_q = calc_Phi_q(joint,q,nBody,nJoint,nCon);
    [A, loc_indx, er] = RowReducedEchelonForm(Phi_q,0);
    m=0;
    for i=(nCon+1):7*nBody
        m=m+1;
        n_prev = n;
        n(m,1) =  loc_indx(1,i);
    end
%     n=4;
    if (n_prev(1,1)~=n(1,1))
        fprintf('Independent variable changed from n=%d to n=%d\n',n_prev(1,1),n(1,1));
        fprintf('v = %f\n',q(n(1,1),1));
        fprintf('v_p = %f\n',q_prime(n(1,1),1));
    end
    % =====================================================================
    % === Adams-Bashforth-Moulton predictor-corrector =====================
    % === first 3 steps with R-K
    % =====================================================================
    if k<4
        [y_out, er] = Num_Int(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,dt);
        m=0;
        for i=1:7*nBody
            if (m<nind && i==n(m+1,1))
                q(i,1) = y_out(m+1,1);
                m=m+1;
            else
            end
        end
        m=0;
        for i=1:7*nBody
            if (m<nind && i==n(m+1,1))
                q_prime(i,1) = y_out(nind+m+1,1);
                m=m+1;
            else
            end
        end
        % update time
        t=t+dt;
        % calculate variables
        [q_out,q_prime_out,q_prime_prime_out,FA,MA,MAT,lambda,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,dt);
        if er ==1
            return;
        end
        % store
        q_store(:,k+1) = q_out;
        q_prime_store(:,k+1) = q_prime_out;
    else % === Adams-Bashforth-Moulton ====================================
        [y_out, er] = Num_Int_2(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,dt,q_store,q_prime_store);
                m=0;
        for i=1:7*nBody
            if (m<nind && i==n(m+1,1))
                q(i,1) = y_out(m+1,1);
                m=m+1;
            else
            end
        end
        m=0;
        for i=1:7*nBody
            if (m<nind && i==n(m+1,1))
                q_prime(i,1) = y_out(nind+m+1,1);
                m=m+1;
            else
            end
        end
        % update time
        t=t+dt;
        % calculate variables
        [q_out,q_prime_out,q_prime_prime_out,FA,MA,MAT,lambda,er] = DIFEQN_3D(q, q_prime, n, joint,nBody,body,nJoint, nCon, t,dt);
        if er ==1
            return;
        end
        % store
        q_store(:,1) = q_store(:,2);
        q_store(:,2) = q_store(:,3);
        q_store(:,3) = q_store(:,4);
        q_store(:,4) = q_out;
        q_prime_store(:,1) =  q_prime_store(:,2);
        q_prime_store(:,2) =  q_prime_store(:,3);
        q_prime_store(:,3) =  q_prime_store(:,4);
        q_prime_store(:,4) =  q_prime_out;
    end
    % =====================================================================
    % ========  UPDATE  ===================================================
    % ======== update Bodies ==============================================
    % =====================================================================
    % redirect for new time step
    q = q_out;
    q_prime = q_prime_out;
    % bodies
    for i=1:nBody
        for j=1:3
            body(i).r(j) = q_out(j+(i-1)*7,1);
            body(i).v(j) = q_prime_out(j+(i-1)*7,1);
        end
        for j=1:4
            body(i).p(j) = q_out(j+(i-1)*7+3,1);
            p_prime{i}(j,1) = q_prime_out(j+(i-1)*7+3,1);
        end
        body(i).A = calc_A(body(i).p');
        body(i).omega = (2*calc_G(body(i).p')*p_prime{i})';
    end
    % =====================================================================
    % ======== EXTRACT FORCES AND MOMENT ON CRANK JOINT
    % === calculate Phi_q
    Phi_q = calc_Phi_q(joint,q_out,nBody,nJoint,nCon);
    % === split Phi_q on Phi_r and Phi_p
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
    hi = (joint(1).q1 - joint(1).p1)';
    fi = (joint(1).r1 - joint(1).p1)';
    gi = -(scew_sym(fi)*hi);
    C_1_REV = [fi,gi,hi];
    lambda_rev = lambda(1:5,1);
    Phi_ri_rev = Phi_r(1:5,1:3);
    % force
    F = (-C_1_REV'*body(1).A'*Phi_ri_rev'*lambda_rev);




    % ======== OUTPUT
    x(k,1) = t;
    x(k,2) = q_out(15,1);
    x(k,3) = q_prime_out(15,1);
    x(k,4) = body(1).omega(1,3);
%     x(k,5) = F(1,1)/100000;

    % ===============================
    % ======== trace and plot tracers
    % ===============================
    tr11 = body(1).r' + body(1).A*t11';
    tr12 = body(1).r' + body(1).A*t12';
    tr21 = body(2).r' + body(2).A*t21';
    tr22 = body(2).r' + body(2).A*t22';
    tr31 = body(3).r' + body(3).A*t31';
    tr32 = body(3).r' + body(3).A*t32';
    tr33 = body(3).r' + body(3).A*t33';
    tr34 = body(3).r' + body(3).A*t34';
    tr41 = body(4).r' + body(4).A*t41';
    tr42 = body(4).r' + body(4).A*t42';
    tr43 = body(4).r' + body(4).A*t43';
    tr44 = body(4).r' + body(4).A*t44';
    % == plot
    crank(1,:) = tr11;
    crank(2,:) = tr12;
    arm(1,:) = tr21;
    arm(2,:) = tr22;
    slider(1,:) = tr31;
    slider(2,:) = tr32;
    slider(3,:) = tr33;
    slider(4,:) = tr34;
    slider(5,:) = tr31;
    ground(1,:) = tr41;
    ground(2,:) = tr42;
    ground(3,:) = tr43;
    ground(4,:) = tr44;
    ground(5,:) = tr41;
    hold on
    clf
    plot3(crank(:,1),crank(:,2),crank(:,3),'o-',arm(:,1),arm(:,2),arm(:,3),'s-',slider(:,1),slider(:,2),slider(:,3),'+-',ground(:,1),ground(:,2),ground(:,3),'+-');
    my_str = ['time = ', num2str(t), 's'];
    title(my_str);
    axis([-3 6 -3 3 -3 3])
    view(0,90)
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    grid on
    drawnow
    pause(0.00001);
    hold off
end




%% =========== PLOTS

%% plot against the book
src1 = 'results/data/x_3.dat';
src2 = 'results/data/x_prime_3_2.dat';
src3 = 'results/data/phi_prime.dat';
src4 = 'results/data/force.dat';
delimiterIn = ' ';
x_book = load(src1);
x_p_book = load(src2);
phi_p_book = load(src3);
f_x_book = load(src4);

figure
hold on
plot(x(:,1),x(:,2),'k',x_book(:,1),x_book(:,2),'r');
axis([0 t 1 6])
title('x_(slider)');
hold off
figure
hold on
plot(x(:,1),x(:,3),'k',x_p_book(:,1),x_p_book(:,2),'r');
axis([0 t -100 100])
title('x''_(slider)');
hold off
figure
hold on
plot(x(:,1),x(:,4),'k',phi_p_book(:,1),phi_p_book(:,2),'r');
title('omega_(crank)');
axis([0 t 0 35])
hold off




