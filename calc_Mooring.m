function [F1,F2,F3,M1,M2,M3,ri1,ri2,ri3,sj1,sj2,sj3,l1,l2,l3,l10,l20,l30] = calc_Mooring(q,q_prime,nBody)
    % strip q and q_prime
    for i=1:nBody
        for j=1:3
            r{i}(j,1) = q(j+(i-1)*7,1);
            r_prime{i}(j,1) = q_prime(j+(i-1)*7,1);
        end
        for j=1:4
            p{i}(j,1)= q(j+(i-1)*7+3,1);
            p_prime{i}(j,1) = q_prime(j+(i-1)*7+3,1);
        end
    end
    % alloc A matrices
    for i=1:nBody
        A{i} = zeros(3,3); 
        A{i} = calc_A(p{i});
    end
    % =====================================================================
    % ============ MOORING ================================================
    % =====================================================================
    % === anchors location in GCS
    xm = 853.87;% [m] radius of anchors from platfrom centreline
    ri1 = [-xm; 0.0; -320.0];
    ri2 = [xm/2; sqrt(3)/2*xm; -320.0];
    ri3 = [xm/2; -sqrt(3)/2*xm; -320.0];
    % === fairleads location in LCS of support
    xf = 5.2;%[m] fairlead radius from centreline of support
    sj1 = [-xf; 0.0; 0.835632828895478e+01];
    sj2 = [xf/2; sqrt(3)/2*xf; 0.835632828895478e+01];
    sj3 = [xf/2; -sqrt(3)/2*xf; 0.835632828895478e+01];
    % === correct for mass centre offset
    s = [1.128758461644648e-02;0;0];
    sj1 = sj1+s;
    sj2 = sj2+s;
    sj3 = sj3+s;
    % === vectors anchor-fairlead
    dij1 = r{2} + A{2}*sj1 - ri1;
    dij2 = r{2} + A{2}*sj2 - ri2;
    dij3 = r{2} + A{2}*sj3 - ri3;
    % === length of dij
    l1 = sqrt(dij1'*dij1);
    l2 = sqrt(dij2'*dij2);
    l3 = sqrt(dij3'*dij3);
    % === INITIAL length of dij
%     l10 = 900;
%     l20 = 900;
%     l30 = 900;
    l10 = sqrt(250*250 + 848.67*848.67);
    l20 = sqrt(250*250 + 848.67*848.67);
    l30 = sqrt(250*250 + 848.67*848.67);
%     l30 = 9.184484836615770e+02;
    % === omega support
    omegaj = (2*calc_G(p{2})*p_prime{2});
    % === time derivative of length l
    l1d = 1/l1*dij1'*(r_prime{2}-A{2}*scew_sym(sj1)*omegaj);
    l2d = 1/l2*dij2'*(r_prime{2}-A{2}*scew_sym(sj2)*omegaj);
    l3d = 1/l3*dij3'*(r_prime{2}-A{2}*scew_sym(sj3)*omegaj);
    % === MAGNITUDE OF FORCE
    % spring and damper constants
    k1 = 384243;
    k2 = 384243;
    k3 = 384243;
    b1 = 384243;
    b2 = 384243;
    b3 = 384243;
%     b1=0.0;
%     b2=0.0;
%     b3=0.0;
    % === REJECT TINY VALUES
    
    % force
    f1m = k1*(l1 - l10) + b1*l1d;
    f2m = k2*(l2 - l20) + b2*l2d;
    f3m = k3*(l3 - l30) + b3*l3d;
    % ===== TURN-OFF PUSHING
%     if l1<l10
%         f1m = 0.0;
%     end
%     if l2<l20
%         f2m = 0.0;
%     end
%     if l3<l30
%         f3m = 0.0;
%     end    
    % === FORCES
    F1 = -f1m/l1*dij1;
    F2 = -f2m/l2*dij2;
    F3 = -f3m/l3*dij3;
    % === MOMENTS
    M1 = -f1m/l1*(scew_sym(sj1)*A{2}'*dij1);
    M2 = -f2m/l2*(scew_sym(sj2)*A{2}'*dij2);
    M3 = -f3m/l3*(scew_sym(sj3)*A{2}'*dij3);


end