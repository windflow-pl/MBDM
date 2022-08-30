function [Fb,Mb,r_b,V,a,b,c,d,s_glob] = calc_Boyancy(q,q_prime,nBody)
    % strip q
    for i=1:nBody
        for j=1:3
            r{i}(j,1) = q(j+(i-1)*7);
        end
        for j=1:4
            p{i}(j,1)= q(j+(i-1)*7+3);
        end
    end
    % strip q_prime
    for i=1:nBody
        for j=1:4
            p_prime{i}(j,1) = q_prime(j+(i-1)*7+3,1);
        end
    end
    % alloc A matrices
    for i=1:nBody
        A{i} = zeros(3,3); 
        A{i} = calc_A(p{i});
    end
    % =====================================================================
    % ============ BOYANCY ================================================
    % =====================================================================
     % === cylinder radius
%     R = 4.62557;%[m]
    R = 4.7;
    % ==== vector normal to the water surface
    vec = A{2}'*[0;0;1];
%     a = vec(1,1);
%     b = vec(2,1);
%     c = vec(3,1);
    a = vec(1,1)/(sqrt(vec(1,1)*vec(1,1)+vec(2,1)*vec(2,1)+vec(3,1)*vec(3,1)));    
    b = vec(2,1)/(sqrt(vec(1,1)*vec(1,1)+vec(2,1)*vec(2,1)+vec(3,1)*vec(3,1)));
    c = vec(3,1)/(sqrt(vec(1,1)*vec(1,1)+vec(2,1)*vec(2,1)+vec(3,1)*vec(3,1)));
    if abs(a)<=(10^(-10)) a=0; end
    if abs(b)<=(10^(-10)) b=0; end
    if abs(c)<=(10^(-10)) c=0; end    
    % ==== vector in LCS from CofG to cylinder axis 
    s = [1.128758461644648e-02;0;0];
%     s = [2;0;0];
%     s=[0.0;0.0;0.0];
    s_glob = A{2}*s;
    % ==== location of the CG shifted to the cylinder axis
    d_tmp = r{2}+A{2}*s;
%     A{2}*s
    d = -d_tmp(3,1)/c;
    % === submerged volume
    V = (d/c+120-7.835632828895478e+01)*pi*R*R;
    % === water density
    rho = 1000;% [kg/m3]
    % =====================================================================
    % === calc boyancy force in GCS
    % =====================================================================
    Fb = [0;0;V*rho*9.81];
    % =====================================================================
    % === calc boyancy centre in LCS
    % =====================================================================
    x = a*pi*R*R*R*R/(4*V*c);
    y = b*pi*R*R*R*R/(4*V*c);
    z = pi*R*R/V*(a*a*R*R/(8*c*c)+b*b*R*R/(8*c*c)+d*d/(2*c*c)-10*10/2);
    r_b = [x;y;z];
    % =====================================================================
    % === calc boyancy moment in LCS
    % =====================================================================
    r_b2 = r_b + s;
    r_b_scew = scew_sym(r_b2);
    Fb_7 = A{2}'*[0;0;V*rho*9.81];
    Mb = r_b_scew*Fb_7;
    
    
end