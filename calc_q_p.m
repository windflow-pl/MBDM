function q_prime = calc_q_p(q,joint,nBody,nJoint,nCon)
     % allocate 
    DPhi = zeros(nCon,7*nBody+1);
    % calc Phi_t
    Phi_t= calc_Phi_t(joint,nJoint,nCon);
    % place -Phi_t as a last column of DPhi
    for i=1:nCon
        DPhi(i,nCon+1) = -1.0*Phi_t(i,1);
    end
    % calc Phi_q
    Phi_q = calc_Phi_q(joint,q,nBody,nJoint,nCon);
        % place Phi_q in DPhi
    for i=1:nCon
        for j=1:7*nBody
            DPhi(i,j) = Phi_q(i,j);
        end
    end
%     DPhi
    % ================ solve
    [A, loc_indx,er] = RowReducedEchelonForm(DPhi,1);
    if er==1           
       return; 
    end
    q_prime = get_solution(A, nCon,7*nBody+1,loc_indx);
end