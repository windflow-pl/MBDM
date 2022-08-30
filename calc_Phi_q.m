function Phi_q = calc_Phi_q(joint,q,nBody,nJoint,nCon)
    % alocate r's and p's
    for i=1:nBody
        r{i} = zeros(3,1);
        p{i} = zeros(4,1);
    end
    % strip q
    for i=1:nBody
        for j=1:3
            r{i}(j,1) = q(j+(i-1)*7);
        end
        for j=1:4
            p{i}(j,1)= q(j+(i-1)*7+3);
        end
    end
    % allocate A matrices
    for i=1:nBody
        A{i} = zeros(3,3); 
        A{i} = calc_A(p{i});
    end
    % allocate G matrices
    for i=1:nBody
        G{i} = zeros(3,3); 
        G{i} = calc_G(p{i});
    end
    % allocate Phi_q
    Phi_q = zeros(nCon,nBody*7);
    % calc Phi_q
    m=0;
    for k=1:nJoint
        if(strcmp(joint(k).type,'R'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            % === find normal vectors to hi
            fi = (joint(k).r1 - joint(k).p1)';
            gi = -(scew_sym(fi)*hi);
            s_jP = joint(k).p2';
            s_iP = joint(k).p1';   
            % ====================== EQN 1-2-3
            Phi_ri = -eye(3);
            Phi_rj = eye(3);
            Phi_pi = 2.0*(A{i}*scew_sym(s_iP))*G{i};
            Phi_pj = -2.0*(A{j}*scew_sym(s_jP))*G{j};
            % update Jakobian
            Phi_q(m+1:m+3,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1:m+3,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+1:m+3,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+1:m+3,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 4
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = -2.0*(hj'*A{j}'*A{i}*scew_sym(fi))*G{i};
            Phi_pj = -2.0*(fi'*A{i}'*A{j}*scew_sym(hj))*G{j};
            % update Jakobian
            Phi_q(m+4,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+4,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+4,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+4,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 5
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = -2.0*(hj'*A{j}'*A{i}*scew_sym(gi))*G{i};
            Phi_pj = -2.0*(gi'*A{i}'*A{j}*scew_sym(hj))*G{j};            
            % update Jakobian
            Phi_q(m+5,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+5,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+5,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+5,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj i j hi hj fi gi s_jP s_iP
            m = m+5;
        end  
        if(strcmp(joint(k).type,'C'))            
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            s_jP = joint(k).p2';
            s_iP = joint(k).p1';
            dij =  r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            % === find normal vectors to hi
            fi = (joint(k).r1 - joint(k).p1)';
            gi = -(scew_sym(fi)*hi);            
            % ====================== EQN 1
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = -2.0*(hj'*A{j}'*A{i}*scew_sym(fi))*G{i};
            Phi_pj = -2.0*(fi'*A{i}'*A{j}*scew_sym(hj))*G{j};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+1,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 2
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = -2.0*(hj'*A{j}'*A{i}*scew_sym(gi))*G{i};
            Phi_pj = -2.0*(gi'*A{i}'*A{j}*scew_sym(hj))*G{j};
            % update Jakobian
            Phi_q(m+2,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+2,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+2,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+2,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 3
            Phi_ri = -fi'*A{i}';
            Phi_rj = fi'*A{i}';
            Phi_pi = 2.0*(fi'*scew_sym(s_iP) - dij'*A{i}*scew_sym(fi))*G{i};
            Phi_pj = 2.0*(-fi'*A{i}'*A{j}*scew_sym(s_jP))*G{j};
            % update Jakobian
            Phi_q(m+3,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+3,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+3,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+3,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 4
            Phi_ri = -gi'*A{i}';
            Phi_rj = gi'*A{i}';
            Phi_pi = 2.0*(gi'*scew_sym(s_iP) - dij'*A{i}*scew_sym(gi))*G{i};
            Phi_pj = 2.0*(-gi'*A{i}'*A{j}*scew_sym(s_jP))*G{j};
            % update Jakobian
            Phi_q(m+4,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+4,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+4,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+4,start_rj+4:start_rj+7) = Phi_pj;            
            m=m+4;
            clear Phi_ri Phi_rj Phi_pi Phi_pj i j hi hj fi gi s_jP s_iP dij
        end
        if(strcmp(joint(k).type,'S'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            % ====================== EQN 1-2-3
            Phi_ri = -eye(3);
            Phi_rj = eye(3);
            Phi_pi = 2.0*(A{i}*scew_sym(s_iP))*G{i};
            Phi_pj = 2.0*(-A{j}*scew_sym(s_jP))*G{j};
            % update Jakobian
            Phi_q(m+1:m+3,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1:m+3,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+1:m+3,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+1:m+3,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj i j s_iP s_jP start_ri start_rj
            m=m+3;
        end        
        if(strcmp(joint(k).type,'RC'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dji =  r{i} + A{i}*s_iP - r{j} - A{j}*s_jP;
            % === find normal vectors to hj
            fj = (joint(k).r2 - joint(k).p2)';
            gj = -(scew_sym(fj)*hj);
            % ====================== EQN 1
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = 2.0*(-hj'*A{j}'*A{i}*scew_sym(hi))*G{i};
            Phi_pj = 2.0*(-hi'*A{i}'*A{j}*scew_sym(hj))*G{j};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+1,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 2
            Phi_ri = fj'*A{j}';
            Phi_rj = -fj'*A{j}';
            Phi_pi = 2.0*(-fj'*A{j}'*A{i}*scew_sym(s_iP))*G{i};
            Phi_pj = 2.0*(fj'*scew_sym(s_jP) - dji'*A{j}*scew_sym(fj))*G{j};
            % update Jakobian
            Phi_q(m+2,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+2,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+2,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+2,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 3
            Phi_ri = gj'*A{j};
            Phi_rj = -gj'*A{j}';
            Phi_pi = 2.0*(-gj'*A{j}'*A{i}*scew_sym(s_iP))*G{i};
            Phi_pj = 2.0*(gj'*scew_sym(s_jP) - dji'*A{j}*scew_sym(gj))*G{j};
            % update Jakobian
            Phi_q(m+3,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+3,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+3,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+3,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj i j start_ri start_rj hi hj s_jP s_iP fj gj dji
            m=m+3;
        end    
        if(strcmp(joint(k).type,'T'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dij =  r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            % === find normal vectors to hi
            fi = (joint(k).r1 - joint(k).p1)';
            gi = -(scew_sym(fi)*hi);
            % === find normal vector to hj
            fj = (joint(k).r2 - joint(k).p2)';
            % ====================== EQN 1
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = 2.0*(-hj'*A{j}'*A{i}*scew_sym(fi))*G{i};
            Phi_pj = 2.0*(-fi'*A{i}'*A{j}*scew_sym(hj))*G{j};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+1,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 2
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = 2.0*(-hj'*A{j}'*A{i}*scew_sym(gi))*G{i};
            Phi_pj = 2.0*(-gi'*A{i}'*A{j}*scew_sym(hj))*G{j};
            % update Jakobian
            Phi_q(m+2,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+2,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+2,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+2,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 3
            Phi_ri = -fi'*A{i}';
            Phi_rj = fi'*A{i}';
            Phi_pi = 2.0*(fi'*scew_sym(s_iP) - dij'*A{i}*scew_sym(fi))*G{i};
            Phi_pj = 2.0*(-fi'*A{i}'*A{j}*scew_sym(s_jP))*G{j};
            % update Jakobian
            Phi_q(m+3,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+3,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+3,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+3,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 4
            Phi_ri = -gi'*A{i}';
            Phi_rj = gi'*A{i}';
            Phi_pi = 2.0*(gi'*scew_sym(s_iP) - dij'*A{i}*scew_sym(gi))*G{i};
            Phi_pj = 2.0*(-gi'*A{i}'*A{j}*scew_sym(s_jP))*G{j};
            % update Jakobian
            Phi_q(m+4,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+4,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+4,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+4,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj
            % ====================== EQN 5
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = 2.0*(-fj'*A{j}'*A{i}*scew_sym(fi))*G{i};
            Phi_pj = 2.0*(-fi'*A{i}'*A{j}*scew_sym(fj))*G{j};
            % update Jakobian
            Phi_q(m+5,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+5,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+5,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+5,start_rj+4:start_rj+7) = Phi_pj;
            clear Phi_ri Phi_rj Phi_pi Phi_pj i j hi hj s_jP s_iP fj fi gi dij
            m=m+5;            
        end
        if(strcmp(joint(k).type,'DI'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dij =  r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            C = joint(k).dist;
            % ====================== EQN 1
            Phi_ri = -2.0*dij';
            Phi_rj = 2.0*dij';
            Phi_pi = 2.0*(2.0*dij'*A{i}*scew_sym(s_iP))*G{i};
            Phi_pj = 2.0*(-2.0*dij'*A{j}*scew_sym(s_jP))*G{j};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+1,start_rj+4:start_rj+7) = Phi_pj;
            clear i j start_ri start_rj s_iP s_jP dij C Phi_ri Phi_rj Phi_pi Phi_pj
            m=m+1;
        end         
        if(strcmp(joint(k).type,'G'))
            i=joint(k).bodies(1,1);
            start_ri = (i-1)*7;
            s_iP = joint(k).p1';
            % ====================== EQN 1-2-3
            Phi_ri = eye(3);
            Phi_pi = -2.0*(A{i}*scew_sym(s_iP))*G{i};
            % update Jakobian
            Phi_q(m+1:m+3,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1:m+3,start_ri+4:start_ri+7) = Phi_pi;
            clear Phi_ri Phi_pi s_iP
            % ====================== EQN 4-5-6
            Phi_ri = zeros(3,3);
            Phi_pi = [zeros(3,1),eye(3)];
            % update Jakobian
            Phi_q(m+4:m+6,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+4:m+6,start_ri+4:start_ri+7) = Phi_pi;
            m=m+6;
            clear i Phi_ri Phi_pi start_ri
        end
        if(strcmp(joint(k).type,'RDRV'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            hi = joint(k).q1';            
            Phi_ri = zeros(1,3);
            Phi_rj = zeros(1,3);
            Phi_pi = -2.0*hi'*G{i};
            Phi_pj = 2.0*hi'*A{i}'*A{j}*G{j};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_rj+1:start_rj+3) = Phi_rj;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi;
            Phi_q(m+1,start_rj+4:start_rj+7) = Phi_pj;     
            clear fj fi hi gi hi Phi_ri Phi_rj Phi_pi Phi_pj i j start_ri start_rj
            m=m+1;
        end
        if(strcmp(joint(k).type,'ABX'))
            i=joint(k).bodies(1,1);
            start_ri = (i-1)*7;
            s_iP = [joint(k).dist;0;0];
            Phi_ri = [1;0;0];
            Phi_pi = -2.0*(A{i}*scew_sym(s_iP))*G{i};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi(1,1:4);
            clear Phi_ri Phi_pi s_iP
            m=m+1;            
        end
        if(strcmp(joint(k).type,'ABY'))
            i=joint(k).bodies(1,1);
            start_ri = (i-1)*7;
            s_iP = [0;joint(k).dist;0];
            Phi_ri = [0;1;0];
            Phi_pi = -2.0*(A{i}*scew_sym(s_iP))*G{i};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi(1,1:4);
            clear Phi_ri Phi_pi s_iP
            m=m+1;            
        end
        if(strcmp(joint(k).type,'ABZ'))
            i=joint(k).bodies(1,1);
            start_ri = (i-1)*7;
            s_iP = [0;0;joint(k).dist];
            Phi_ri = [0;0;1];
            Phi_pi = -2.0*(A{i}*scew_sym(s_iP))*G{i};
            % update Jakobian
            Phi_q(m+1,start_ri+1:start_ri+3) = Phi_ri;
            Phi_q(m+1,start_ri+4:start_ri+7) = Phi_pi(3,1:4);
            clear Phi_ri Phi_pi s_iP 
            m=m+1;            
        end
        if(strcmp(joint(k).type,'ANGFIX'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            e02 = joint(k).q2(2);
            e12 = joint(k).r2(1);
            e22 = joint(k).r2(2);
            e32 = joint(k).r2(3);
            % update Jakobian
            Phi_q(m+1:m+3,start_ri+1:start_ri+3) = zeros(3,3);
            Phi_q(m+1:m+3,start_rj+1:start_rj+3) = zeros(3,3);
            Phi_q(m+1:m+3,start_ri+4:start_ri+7) = [zeros(3,1), eye(3)];
            Phi_q(m+1:m+3,start_rj+4:start_rj+7) = [-e12,-e02,e32,-e22; -e22,-e32,-e02,e12; -e32,e22,-e12,-e02];
            m=m+3;
        end
        
        
    end
    %% add math constraints Jakobian
    for i=1:nBody        
        start_ri = (i-1)*7; 
        Phi_ri = zeros(1,3);
        Phi_pi = 2.0*p{i}';
        % update Jakobian
        for j=1:3
            Phi_q(m+i,start_ri+j) = Phi_ri(j);
        end
        for j=1:4
            Phi_q(m+i,start_ri+j+3) = Phi_pi(j);
        end
        clear Phi_ri Phi_pi
    end
end