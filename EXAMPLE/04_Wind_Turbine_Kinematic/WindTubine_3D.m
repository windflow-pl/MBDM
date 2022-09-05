clc
clear all

%% INSTALL MBDM FUNCTIONS
addpath('../../.')

%% input
nBody = 7;
nJoint = 13; % include driving

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
angle = [0, 0, 0];%euler angles 1-2-3 notation
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
body(1).name = 'blade1';
body(1).r = [-5.0,0.0,111.975];%[m]
body(1).p = [e0,e1,e2,e3];
body(1).v = [0,0,0];%[m/s]
body(1).omega = [-121/300*pi,0,0];%[rad/s]
body(1).m = 17740;%[kg]
body(1).I = [4338984,0,0; 0,4338984,0; 0,0,1084746];% kg*m^2
body(1).A = calc_A(body(1).p');
% ===============================
% === inital orientation
angle = [4/3*pi, 0.0, 0.0];
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
body(2).name = 'blade2';
body(2).r = [-5,21.975*sqrt(3)/2,79.0125];
body(2).p = [e0,e1,e2,e3];
% body(2).p = [-0.5,sqrt(3)/2,0,0];
body(2).v = [0,0,0];
body(2).omega = [-121/300*pi,0,0];
body(2).m = 17740;% kg
body(2).I = [4338984,0,0; 0,4338984,0; 0,0,1084746];% kg*m^2
body(2).A = calc_A(body(2).p');
% ===============================
% === inital orientation
angle = [2/3*pi, 0.0, 0.0];
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
body(3).name = 'blade3';
body(3).r = [-5,-21.975*sqrt(3)/2,79.0125];
body(3).p = [e0,e1,e2,e3];
% body(2).p = [0.5,sqrt(3)/2,0,0];
body(3).v = [0,0,0];
body(3).omega = [-121/300*pi,0,0];
body(3).m = 17740;% kg
body(3).I = [4338984,0,0; 0,4338984,0; 0,0,1084746];% kg*m^2
body(3).A = calc_A(body(3).p');
% ===============================
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
body(4).name = 'tower';
body(4).r = [0.0,0.0,43.4];
body(4).p = [e0,e1,e2,e3];
body(4).v = [0,0,0];
body(4).omega = [0,0,0];
body(4).m = 249718;% kg
body(4).I = [31987367,0,0; 31987367,0,0; 0,0,1318823];% kg*m^2
body(4).A = calc_A(body(4).p');
% ===============================
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
body(5).name = 'nacelle';
body(5).r = [1.9,0.0,90.0];
body(5).p = [e0,e1,e2,e3];
body(5).v = [0,0,0];
body(5).omega = [0,0,0];
body(5).m = 240000;% kg
body(5).I = [1741490,0,0; 0,1741490,0; 0,0,1741490];% kg*m^2
body(5).A = calc_A(body(5).p');
% ===============================
% === inital orientation
angle = [0, 0.0, 0.0];
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
body(6).name = 'spinner';
body(6).r = [-5,0.0,90.0];
body(6).p = [e0,e1,e2,e3];
body(6).v = [0,0,0];
body(6).omega = [-121/300*pi,0,0];
body(6).m = 56780;% kg
body(6).I = [115926,0,0; 0,75171,0; 0,0,75171];% kg*m^2
body(6).A = calc_A(body(6).p');
% ===============================
% === inital orientation
angle = [0, 0.0, 0.0];
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
body(7).name = 'support';
body(7).r = [0.0,0.0,-89.9155];
body(7).p = [e0,e1,e2,e3];
body(7).v = [0,0,0];
body(7).omega = [0,0,0];
body(7).m = 7466330;% kg
body(7).I = [4229230000,0,0; 0,4229230000,0; 0,0,164230000];% kg*m^2
body(7).A = calc_A(body(6).p');
% ====== clear
clear angle A e0 e1 e2 e3
% ===============================
% ###############################
% ======== define joints ========
% ###############################
% ===============================
%%%%%%%%%%% BLADE 1 %%%%%%%%%%%
% % revolute joint
joint(1).type = 'R';
joint(1).bodies = [1,6];%IDs of body 1 and body 2 (4 tower 5 nacelle)
joint(1).p1 = [0,0,-20.475];% joint location in LCS of body 1
joint(1).q1 = [0,1,-20.475];% point on the axis of rotation in LCS 1
joint(1).r1 = [0,0,-19.475];% point located on the normal to the axis of rotation in LCS 1
joint(1).p2 = [0,0,1.5];% joint location in LCS of body 2
joint(1).q2 = [0,1,1.5];% point on the axis of rotation in LCS 2
% === % RDRV for nacelle
joint(2).type = 'RDRV';% rotational driving
joint(2).bodies = [6,1];%ID 2 is ID of the body to be driven relative to body 1
joint(2).p2 = [0,0,1];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
% angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
joint(2).p1 = [0,0,1];% versor in LCS 1 that lies in plane of rotation (see above)
joint(2).q1 = [0,1,0];% versor in LCS 1 indicating the axis of rotation
joint(2).r1 = [0,0,0];% inital angle in the first entry (from p1 to p2) 0<angle<2*pi
joint(2).dist = 0;% rotational velocity (positive counter-clockwise)
%%%%%%%%%%% BLADE 2 %%%%%%%%%%%
% % revolute joint
joint(3).type = 'R';
joint(3).bodies = [2,6];%IDs of body 1 and body 2 (4 tower 5 nacelle)
joint(3).p1 = [0,0,-20.475];% joint location in LCS of body 1
joint(3).q1 = [1,0,-20.475];% point on the axis of rotation in LCS 1
joint(3).r1 = [0,0,-19.475];% point located on the normal to the axis of rotation in LCS 1
joint(3).p2 = [0,1.5*sqrt(3)/2,-0.75];% joint location in LCS of body 2
joint(3).q2 = [1,1.5*sqrt(3)/2,-0.75];% point on the axis of rotation in LCS 2
% === % RDRV for nacelle
joint(4).type = 'RDRV';% rotational driving
joint(4).bodies = [6,2];%ID 2 is ID of the body to be driven relative to body 1
joint(4).p2 = [0,0,1];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
% angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
joint(4).p1 = [0,0,1];% versor in LCS 1 that lies in plane of rotation (see above)
joint(4).q1 = [1,0,0];% versor in LCS 1 indicating the axis of rotation
joint(4).r1 = [4/3*pi,0,0];% inital angle in the first entry (from p1 to p2) 0<angle<2*pi
joint(4).dist = 0;% rotational velocity (positive counter-clockwise)
%%%%%%%%%%% BLADE 3 %%%%%%%%%%%
% % revolute joint
joint(5).type = 'R';
joint(5).bodies = [3,6];%IDs of body 1 and body 2 (4 tower 5 nacelle)
joint(5).p1 = [0,0,-20.475];% joint location in LCS of body 1
joint(5).q1 = [1,0,-20.475];% point on the axis of rotation in LCS 1
joint(5).r1 = [0,0,-19.475];% point located on the normal to the axis of rotation in LCS 1
joint(5).p2 = [0,-1.5*sqrt(3)/2,-0.75];% joint location in LCS of body 2
joint(5).q2 = [1,-1.5*sqrt(3)/2,-0.75];% point on the axis of rotation in LCS 2
% === % RDRV for nacelle
joint(6).type = 'RDRV';% rotational driving
joint(6).bodies = [6,3];%ID 2 is ID of the body to be driven relative to body 1
joint(6).p2 = [0,0,1];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
% angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
joint(6).p1 = [0,0,1];% versor in LCS 1 that lies in plane of rotation (see above)
joint(6).q1 = [1,0,0];% versor in LCS 1 indicating the axis of rotation
joint(6).r1 = [2/3*pi,0,0];% inital angle in the first entry (from p1 to p2) 0<angle<2*pi
joint(6).dist = 0;% rotational velocity (positive counter-clockwise)
%%%%%%%%%%% NACELLE-TOWER %%%%%%%%%%%
% % revolute joint
joint(7).type = 'R';
joint(7).bodies = [4,5];%IDs of body 1 and body 2 (4 tower 5 nacelle)
joint(7).p1 = [0,0,46.6];% joint location in LCS of body 1
joint(7).q1 = [0,0,47.6];% point on the axis of rotation in LCS 1
joint(7).r1 = [0,1,46.6];% point located on the normal to the axis of rotation in LCS 1
joint(7).p2 = [-1.9,0.0,0];% joint location in LCS of body 2
joint(7).q2 = [-1.9,0.0,1];% point on the axis of rotation in LCS 2
% === % RDRV for nacelle
joint(8).type = 'RDRV';% rotational driving
joint(8).bodies = [5,4];%ID 2 is ID of the body to be driven relative to body 1
joint(8).p2 = [1,0,0];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
% angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
joint(8).p1 = [1,0,0];% versor in LCS 1 that lies in plane of rotation (see above)
joint(8).q1 = [0,0,1];% versor in LCS 1 indicating the axis of rotation
joint(8).r1 = [0,0,0];% inital angle in the first entry (from p1 to p2) 0<angle<2*pi
joint(8).dist = 0;% rotational velocity (positive counter-clockwise)
%%%%%%%%%%% TOWER-SUPPORT %%%%%%%%%%%
% === % ANGFIX (angle fix) for tower
joint(9).type = 'ANGFIX';
joint(9).bodies = [4,7];%IDs of body 1 and body 2 (1 will be fixed to 2) % 4 tower 7 support
angle = [0.0, 0.0, 0.0]; % 1-2-3 notation for body 1 orientation relative to body 2
A = calc_A_angle(angle);
e0 = sqrt((trace(A)+1)/4);
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
joint(9).p1 = [body(4).p(2),body(4).p(3),body(4).p(4)];%
joint(9).p2 = [body(7).p(2),body(7).p(3),body(7).p(4)];%
joint(9).q2 = [body(7).p(1),e0,0];% first entry e0 body 2, second e0 of additional rotation
joint(9).r2 = [e1,e2,e3];% Euler parameres orienting body 1 elative to 2 [e1,e2,e3];
% === % spherical joint fixing distance for nacelle
joint(10).type = 'S';
joint(10).bodies = [4,7];%IDs of body 1 and body 2
joint(10).p1 = [0,0,-33.4];% joint location in LCS of body 1
joint(10).p2 = [0,0,99.9155];% joint location in LCS of body 2
% % revolute joint
joint(11).type = 'R';
joint(11).bodies = [5,6];%IDs of body 1 and body 2 (5 nacelle, 6 spinner)
joint(11).p1 = [-6.3375,0,0];% joint location in LCS of body 1
joint(11).q1 = [-5.3375,0,0];% point on the axis of rotation in LCS 1
joint(11).r1 = [-6.3375,1,0];% point located on the normal to the axis of rotation in LCS 1
joint(11).p2 = [0.5625,0.0,0.0];% joint location in LCS of body 2
joint(11).q2 = [1.5625,0.0,0.0];% point on the axis of rotation in LCS 2
% % ground
joint(12).type = 'G';
joint(12).bodies = [7,0];%ID of the body to be fixed to the ground
joint(12).p1 = [0,0,-89.9155];% desired location of LCS in GCS
joint(12).q1 = [0,0,0];% desired orientation in terms of Euler parameters [p2,p3,p4]
% === % RDRV for spinner
joint(13).type = 'RDRV';% rotational driving
joint(13).bodies = [6,5];%ID 2 is ID of the body to be driven relative to body 1
joint(13).p2 = [0,1,0];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
% angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
joint(13).p1 = [0,1,0];% versor in LCS 1 that lies in plane of rotation (see above)
joint(13).q1 = [1,0,0];% versor in LCS 1 indicating the axis of rotation
joint(13).r1 = [0,0,0];% inital angle in the first entry (from p1 to p2) 0<angle<2*pi
joint(13).dist = -8.836*2*pi/60;% rotational velocity (positive counter-clockwise)

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
    if(strcmp(joint(i).type,'ANGFIX')) nConstraints =nConstraints+3; end
end
disp(sprintf('Number of bodies = %d, Joints + Driving = %d',nBody,nJoint));
disp(sprintf('Number of constraint equations = %d, DOFs = %d',nConstraints,DOFs-nConstraints));
% == number of constraints with quaternions
nCon = nConstraints+nBody;

%% =========    KINEMATICS     ==========
% =======================================
er =0;
t = 0.000;
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

% =======================================
% ======= define time step
% =======================================
dt = 0.1;
t_final = 20;

% =======================================
% ======= define tracers
% =======================================
t11 = [0.0, 0.0, -20.475];
t12 = [0.0, 0.0, 41.025];
t21 = [0.0, 0.0, -20.475];
t22 = [0.0, 0.0, 41.025];
t31 = [0.0, 0.0, -20.475];
t32 = [0.0, 0.0, 41.025];
t41 = [0,0,-33.4];% tower
t42 = [0,0,46.6];
t51 = [0,0,0];% nacelle
t52 = [-6.3375,0,0];
t61 = [0,0,0];% spinner
t62 = [0.5625,0,0];
t71 = [0,0,80];% support
t72 = [0,0,99.9155];

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
    [q_out,Phi_q,er] = N_R_kinematic(q, joint,nBody,nJoint, nCon, t);
%     rank(Phi_q)
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
            p_prime_prime{i}(j,1) = q_prime_prime(j+(i-1)*7+3,1);
        end
        body(i).A = calc_A(body(i).p');
        body(i).omega = (2*calc_G(body(i).p')*p_prime{i})';
    end
    % === create output of variables
    omega_prime = (2*calc_G(body(1).p')*p_prime_prime{1})';
    x(k,1) = t;
    x(k,2) = body(1).omega(1);
    x(k,3) = body(2).omega(1);
    x(k,4) = body(3).omega(1);
    % =====================================================================
    % ======== END OF KINEMATICS ==========================================
    % =====================================================================

    % ===============================
    % ======== trace and plot tracers
    % ===============================
    Tr11 = body(1).r' + body(1).A*t11';
    Tr12 = body(1).r' + body(1).A*t12';
    Tr21 = body(2).r' + body(2).A*t21';
    Tr22 = body(2).r' + body(2).A*t22';
    Tr31 = body(3).r' + body(3).A*t31';
    Tr32 = body(3).r' + body(3).A*t32';
    Tr41 = body(4).r' + body(4).A*t41';
    Tr42 = body(4).r' + body(4).A*t42';
    Tr51 = body(5).r' + body(5).A*t51';
    Tr52 = body(5).r' + body(5).A*t52';
    Tr61 = body(6).r' + body(6).A*t61';
    Tr62 = body(6).r' + body(6).A*t62';
    Tr71 = body(7).r' + body(7).A*t71';
    Tr72 = body(7).r' + body(7).A*t72';

    blade_1(1,:) = Tr11;
    blade_1(2,:) = Tr12;
    blade_2(1,:) = Tr21;
    blade_2(2,:) = Tr22;
    blade_3(1,:) = Tr31;
    blade_3(2,:) = Tr32;
    tower(1,:) = Tr41;
    tower(2,:) = Tr42;
    nacelle(1,:) = Tr51;
    nacelle(2,:) = Tr52;
    spinner(1,:) = Tr61;
    spinner(2,:) = Tr62;
    support(1,:) = Tr71;
    support(2,:) = Tr72;


    h1 = figure(1);
    set(h1,'Position',[10 10 500 500])
    set(h1, 'DefaultLineLineWidth', 2);
    clf;
    subplot (2,2,1, "align");
    plot3(blade_1(:,1),blade_1(:,2),blade_1(:,3),'o-',blade_2(:,1),blade_2(:,2),blade_2(:,3),'o-',blade_3(:,1),blade_3(:,2),blade_3(:,3),'o-',tower(:,1),tower(:,2),tower(:,3),'s-',nacelle(:,1),nacelle(:,2),nacelle(:,3),'+-',spinner(:,1),spinner(:,2),spinner(:,3),'+-',support(:,1),support(:,2),support(:,3),'s-');
    axis([-150 150 -150 150 -10 170])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    view(-90,0)
    title(['t=', num2str(t,'%1.5f\n')])
    daspect([1 1 1])
    grid on

    subplot (2,2,2, "align");
    plot3(blade_1(:,1),blade_1(:,2),blade_1(:,3),'o-',blade_2(:,1),blade_2(:,2),blade_2(:,3),'o-',blade_3(:,1),blade_3(:,2),blade_3(:,3),'o-',tower(:,1),tower(:,2),tower(:,3),'s-',nacelle(:,1),nacelle(:,2),nacelle(:,3),'+-',spinner(:,1),spinner(:,2),spinner(:,3),'+-',support(:,1),support(:,2),support(:,3),'s-');
    axis([-150 150 -150 150 -10 170])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    view(0,0)
    title(['t=', num2str(t,'%1.5f\n')])
    daspect([1 1 1])
    grid on

    subplot (2,2,[3 4], "align");
    plot3(blade_1(:,1),blade_1(:,2),blade_1(:,3),'o-',blade_2(:,1),blade_2(:,2),blade_2(:,3),'o-',blade_3(:,1),blade_3(:,2),blade_3(:,3),'o-',tower(:,1),tower(:,2),tower(:,3),'s-',nacelle(:,1),nacelle(:,2),nacelle(:,3),'+-',spinner(:,1),spinner(:,2),spinner(:,3),'+-',support(:,1),support(:,2),support(:,3),'s-');
    axis([-150 150 -150 150 -10 170])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    view(-50,20)
    title(['t=', num2str(t,'%1.5f\n')])
    daspect([1 1 1])
    grid on

    %save image
    figname = sprintf('./results/img/figure_%04d.png',k);
    print ("-r100", figname);


end


