function q_prime_prime = calc_q_pp(q,q_prime,joint,nBody,nJoint,nCon)
    % allocate 
    DPhi = zeros(nCon,nCon+1);
    % calc Gamma
    gamma = calc_gamma(joint,q,q_prime,nBody,nJoint,nCon);    
    % place Gamma as a last column of DPhi
    for i=1:nCon
        DPhi(i,7*nBody+1) = gamma(i,1);
    end
    % calc Phi_q
    Phi_q = calc_Phi_q(joint,q,nBody,nJoint,nCon);
    % place Phi_q in DPhi
    for i=1:nCon
        for j=1:7*nBody
            DPhi(i,j) = Phi_q(i,j);
        end
    end
    % ================ solve
    [A, loc_indx,er] = RowReducedEchelonForm(DPhi,1);
    if er==1           
       return; 
    end
    q_prime_prime = get_solution(A, nCon,7*nBody+1,loc_indx);
    
end