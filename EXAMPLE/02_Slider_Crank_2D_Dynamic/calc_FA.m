function FA = calc_FA(q,q_prime,nBody,body)
    % allocate
    FA = zeros(3*nBody,1);
    for i=1:nBody
        m=(i-1)*3;
        FA(m+1:m+3,1) = [body(i).m*9.81;0;0];        
    end
    x3 = q(15,1);
    v3 = q_prime(15,1);
    
%     x3_FA = x3
%     v3_FA = v3
    F = 0.0;
    if(v3>0)
        if(1.5<=x3 && x3<=5.0)
            F = -282857/(6-x3) + 62857;
        elseif(5.0<x3 && x3<=5.6)
            F = -110000*(1-sin(2*pi*(x3-5.25)));
        end
    end    
    if(v3<=0)
        F = 0.0;
    end
%     F
    FA(7,1) = FA(7,1) + F;
end