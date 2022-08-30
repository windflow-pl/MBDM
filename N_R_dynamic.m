function [q_out,DPhi,Phi_u,Phi_v,er] =  N_R_dynamic(q,q_prime,n,joint,nBody,nJoint, nCon, t,dt)
% [q_out,er] = N_R_dynamic(q,n,joint,nBody,nJoint, nCon, t)
    tolq=1.0e-5;
    tolf=1.0e-5;
    steps = 200;
    % number of independent variables
    nind = 7*nBody-nCon;
    % ALLOCATE
    DPhi = zeros(7*nBody,7*nBody+1);
    % redirect input to output before going for iteration
    er = 0;
    dspld=0;
    q_out = q;
    % ==== correct q with velocities to get closer to solution
%     m=0;
%     for i=1:7*nBody      
%         if (m<nind && i==n(m+1,1))             
%         else
%             q_out(i,1) = q_out(i,1) + q_prime(i,1)*dt;
%         end
%     end
    % allocate Phi_u
    Phi_u = zeros(nCon,nCon);
    % allocate Phi_v
    Phi_v = zeros(nCon,7*nBody-nCon);
    % allocate u and v
    u = zeros(nCon,1);
    v = zeros(7*nBody-nCon,1);    
    for k=1:steps 
        % calc Phi
        [Phi, er] = calc_Phi(joint,q_out,t,nBody,nJoint,nCon);
        if er==1
           fprintf('ERROR: angle has imaginary part in N-R\n');
           return; 
        end
        % check if Phi[i]<tolf
        j=0;
        for i=1:nCon
            if (abs(Phi(i,1))>=tolf) j=1; end
        end
        if (j==0) 
%             sprintf('N-R reached conv in tolf\n')
            break; 
        end
        % place [-Phi] as a last column of DPhi
        for i=1:nCon           
            DPhi(i,7*nBody+1) = -1.0*Phi(i,1);
        end
        % calc Phi_q
        Phi_q = calc_Phi_q(joint,q_out,nBody,nJoint,nCon);
        % extract Phi_u and Phi_v
        m=0;
        l=0;
        for i=1:7*nBody
            if (m<nind && i==n(m+1,1))
%                 n(m+1,1)
                Phi_v(:,m+1) = Phi_q(:,i);
                m=m+1;
            else
                l=l+1;
                Phi_u(:,l) = Phi_q(:,i); 
            end
        end
        % extract u and v
        m=0;
        l=0;
        for i=1:7*nBody
            if (m<nind && i==n(m+1,1))
                v(m+1,1) = q_out(i,1);
                m=m+1;
            else
                l=l+1;
                u(l,1) = q_out(i,1); 
            end
        end
        % place Phi_u and Phi_v in DPhi
        for i=1:nCon
            for j=1:nCon
                DPhi(i,j) = Phi_u(i,j);
            end
        end
        for i=1:nCon
            for j=1:(7*nBody-nCon)
                DPhi(i,j+nCon) = Phi_v(i,j);
            end
        end
        % fill last rows of DPhi with 0 (already filled) and I matrix
        for i=(nCon+1):7*nBody
            for j=(nCon+1):7*nBody
                if i==j
                    DPhi(i,j) = 1.0;
                else
                    DPhi(i,j) = 0.0;
                end
            end
        end
        % ================ solve
        [A, loc_indx, er] = RowReducedEchelonForm(DPhi,1);
        if er==1           
           return; 
        end
        dsol = get_solution(A, 7*nBody,7*nBody+1,loc_indx); 
        % check if dsol<tolq
        j=0;
        for i=1:nCon
            if(abs(dsol(i)) >=tolq) j=1; end  
        end        
        if (j==0) 
%             sprintf('N-R reached conv in tolq\n') 
            break; 
        end        
        % UPDATE
        dep=0;
        ind=0;
        for i=1:7*nBody
            if (ind<nind && i==n(ind+1,1))
                ind=ind+1;
                q_out(i,1) = q_out(i,1) + dsol(nCon+ind,1);                
            else
                dep=dep+1;
                q_out(i,1) = q_out(i,1) + dsol(dep,1);
            end
        end
        
%         if k>10 && dspld == 0
%             fprintf('N-R went over 10 stpes...\n');
%             dspld=1;
% %             dsol
% %             Phi = calc_Phi(joint,q_out,t,nBody,nJoint,nCon);
% %             Phi
%         end
        
        
        if k==steps
            er = 1;
            disp(sprintf('N-R didnt converge in %d steps. Try to increase ''steps'' or check Jakobian.\n',k));
            return;
        end  
    end    
end