function [Thrust,x] = calc_Thrust(q,nBody,blade_ID,hub_ID,nacelle_ID)
    a0 = 5.471261419731895e+02;
    a1 = -9.224386003387494e-01;
    b1 = 7.481387825666774e+00;
    a2 = 5.687687399803831e+00;
    b2 = 2.838357616988216e-02;
    a3 = -2.206893100889292e+00;
    b3 = -3.708819643504395e+00;
    a4 = -1.725334280068869e+00;
    b4 = 2.333835354164476e+00;
    a5 = 1.677366004171450e+00;
    b5 = 6.277033100479327e-01;
    w = 2.445983032762768e-02;
    % === initiate vesrors (as in constrains)
    p2 = [0,0,1];% versor in LCS 2 that lies in plane of rotation (is a reference to calculate angle!)
    % angle is calculated as angle >>between p1 and p2<< (from p1 to p2). If p1=p2 => angle=0
    p1 = [0,0,1];% versor in LCS 1 that lies in plane of rotation (see above)
    q1 = [1,0,0];% versor in LCS 1 indicating the axis of rotation        
    % === strip q
    for i=1:nBody
        for j=1:3
            r{i}(j,1) = q(j+(i-1)*7);
        end
        for j=1:4
            p{i}(j,1)= q(j+(i-1)*7+3);
        end
    end
    for i=1:nBody
        A{i} = zeros(3,3); 
        A{i} = calc_A(p{i});
    end
    % === specify indeces
    j = hub_ID;
    i = nacelle_ID;
    % === find angle of the blade-hub
    fj = p2';
    fi = p1';
    hi = q1';
    gi = -(scew_sym(fi)*hi);
    c = fi'*A{i}'*A{j}*fj;
    s = gi'*A{i}'*A{j}*fj; 
    % ==== fix cos and sin if too big or too small
    % quaternions are normalized up to accuracy of N-R
    % and transformation matrices A can scale versors
    % if |p|>1 which makes cos and sin too big or too small
    if c>1
        c=1.0;
        s=0.0;
    end
    if s>1
        s=1.0;
        c=0.0;
    end
    if c<-1
        c=-1.0;
        s = 0.0;
    end
    if s<-1
        s=-1.0;
        c = 0.0;
    end
    % == find angle of hub-nacelle
    if (s>=0 && c>=0)
        ang = asin(s);
    elseif (s>=0 && c<0)
        ang = pi - asin(s);
    elseif (s<0 && c<0)
        ang = pi - asin(s);
    elseif (s<0 && c>=0)
        ang = 2*pi + asin(s);
    end
    % === add angle hub-blade
    if blade_ID==1
        x = ang;
    end
    if blade_ID==2
        x = ang+4/3*pi;
    end
    if blade_ID==3
        x = ang+2/3*pi;
    end
    % === Thrust [N]
    Thrust  =  156.934283047*(a0+a1*cos(x*w)+b1*sin(x*w)+a2*cos(2*x*w)+b2*sin(2*x*w)+a3*cos(3*x*w)+b3*sin(3*x*w)+a4*cos(4*x*w)+b4*sin(4*x*w)+a5*cos(5*x*w)+b5*sin(5*x*w));
end