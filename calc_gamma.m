function gamma = calc_gamma(joint,q,q_prime,nBody,nJoint,nCon)
    % == alocate r's and p's
    for i=1:nBody
        r{i} = zeros(3,1);
        p{i} = zeros(4,1);
    end
    % == strip q
    for i=1:nBody
        for j=1:3
            r{i}(j,1) = q(j+(i-1)*7);
        end
        for j=1:4
            p{i}(j,1)= q(j+(i-1)*7+3);
        end
    end
    % == alocate r_prime and p_prime
    for i=1:nBody
        r_prime{i} = zeros(3,1);
        p_prime{i} = zeros(4,1);
    end
    % == strip q_prime
    for i=1:nBody
        for j=1:3
            r_prime{i}(j,1) = q_prime(j+(i-1)*7);
        end
        for j=1:4
            p_prime{i}(j,1)= q_prime(j+(i-1)*7+3);
        end
    end
    % == allocate A matrices
    for i=1:nBody
        A{i} = zeros(3,3); 
        A{i} = calc_A(p{i});
    end
    % == allocate G matrices
    for i=1:nBody
        G{i} = zeros(3,3); 
        G{i} = calc_G(p{i});
    end
    % == alocate omegas (defined in LCS)
    for i=1:nBody
        omega{i} = zeros(3,1);
    end
    % == calculate omegas
    for i=1:nBody
        omega{i} = 2.0*G{i}*p_prime{i}; 
    end
    % == allocate gamma
    gamma = zeros(nCon,1);
    % ==================================
    % === main loop
    % ==================================
    m=0;
    for k=1:nJoint
        if(strcmp(joint(k).type,'R'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            % === find normal vectors to hi
            fi = (joint(k).r1 - joint(k).p1)';
            gi = -(scew_sym(fi)*hi);
            s_jP = joint(k).p2';
            s_iP = joint(k).p1'; 
            ss_omega_i = scew_sym(omega{i});
            ss_omega_j = scew_sym(omega{j});
            % ====================== EQN 1-2-3
            gamma(m+1:m+3,1) = A{i}*ss_omega_i*ss_omega_i*s_iP - A{j}*ss_omega_j*ss_omega_j*s_jP;
            % ====================== EQN 4
            gamma(m+4,1) = -hj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*fi + 2*omega{j}'*scew_sym(hj)*A{j}'*A{i}*scew_sym(fi)*omega{i};
            % ====================== EQN 5
            gamma(m+5,1) = -hj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*gi + 2*omega{j}'*scew_sym(hj)*A{j}'*A{i}*scew_sym(gi)*omega{i};
            m=m+5;
            clear i j hi hj fi gi s_jP s_iP ss_omega_i ss_omega_j
        end
        if(strcmp(joint(k).type,'C'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            start_ri = (i-1)*7;
            start_rj = (j-1)*7;
            ss_omega_i = scew_sym(omega{i});
            ss_omega_j = scew_sym(omega{j});
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            s_jP = joint(k).p2';
            s_iP = joint(k).p1';
            dij =  r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            % === find normal vectors to hi
            fi = (joint(k).r1 - joint(k).p1)';
            gi = -(scew_sym(fi)*hi); 
            % ====================== EQN 1
            gamma(m+1,1) = -hj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*fi + 2*omega{j}'*scew_sym(hj)*A{j}'*A{i}*scew_sym(fi)*omega{i};
            % ====================== EQN 2
            gamma(m+2,1) = -hj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*gi + 2*omega{j}'*scew_sym(hj)*A{j}'*A{i}*scew_sym(gi)*omega{i};
            % ====================== EQN 3
            tmp1 = 2*omega{i}'*scew_sym(fi)*A{i}'*(r_prime{i}-r_prime{j}) + 2*s_jP'*ss_omega_j*A{j}'*A{i}*fi;
            tmp2 = -s_iP'*ss_omega_i*ss_omega_i*fi - s_jP'*ss_omega_j*ss_omega_j*A{j}'*A{i}*fi - dij'*A{i}*ss_omega_i*ss_omega_i*fi;
            gamma(m+3,1) = tmp1 + tmp2;
            clear tmp1 tmp2
            % ====================== EQN 4
            tmp1 = 2*omega{i}'*scew_sym(gi)*A{i}'*(r_prime{i}-r_prime{j}) + 2*s_jP'*ss_omega_j*A{j}'*A{i}*gi;
            tmp2 = -s_iP'*ss_omega_i*ss_omega_i*gi - s_jP'*ss_omega_j*ss_omega_j*A{j}'*A{i}*gi - dij'*A{i}*ss_omega_i*ss_omega_i*gi;
            gamma(m+4,1) = tmp1 + tmp2;
            clear tmp1 tmp2
            m=m+4;
            clear i j hi hj s_iP s_jP dij ss_omega_i ss_omega_j fi gi fj
        end
        if(strcmp(joint(k).type,'S'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            ss_omega_i = scew_sym(omega{i});
            ss_omega_j = scew_sym(omega{j});
            % ====================== EQN 1-2-3
            gamma(m+1:m+3,1) = A{i}*ss_omega_i*ss_omega_i*s_iP - A{j}*ss_omega_j*ss_omega_j*s_jP;
            m=m+3;
            clear i j s_iP s_jP ss_omega_i ss_omega_j
        end
        if(strcmp(joint(k).type,'RC'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dji =  r{i} + A{i}*s_iP - r{j} - A{j}*s_jP;
            ss_omega_i = scew_sym(omega{i});
            ss_omega_j = scew_sym(omega{j});
            % === find normal vectors to hj
            fj = (joint(k).r2 - joint(k).p2)';
            gj = -(scew_sym(fj)*hj);
            % ====================== EQN 1
            gamma(m+1,1) = -hj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*hi + 2*omega{j}'*scew_sym(hj)*A{j}'*A{i}*scew_sym(hi)*omega{i};
            % ====================== EQN 2            
            tmp1 = 2*omega{j}'*scew_sym(fj)*A{j}*(r_prime{j}-r_prime{i}) + 2*s_iP'*ss_omega_i*A{i}'*A{j}*ss_omega_j*fj;
            tmp2 = -s_iP'*ss_omega_j*ss_omega_j*fj - s_iP'*ss_omega_i*ss_omega_i*A{i}'*A{j}*fj - dji'*A{j}*ss_omega_j*ss_omega_j*fj;
            gamma(m+2,1) = tmp1 + tmp2;
            clear tmp1 tmp2
            % ====================== EQN 3
            tmp1 = 2*omega{j}'*scew_sym(gj)*A{j}*(r_prime{j}-r_prime{i}) + 2*s_iP'*ss_omega_i*A{i}'*A{j}*ss_omega_j*gj;
            tmp2 = -s_iP'*ss_omega_j*ss_omega_j*gj - s_iP'*ss_omega_i*ss_omega_i*A{i}'*A{j}*gj - dji'*A{j}*ss_omega_j*ss_omega_j*gj;
            gamma(m+3,1) = tmp1 + tmp2;
            m=m+3;
            clear i j hi hj s_iP s_jP dji ss_omega_i ss_omega_j fj gj tmp1 tmp2
        end
        if(strcmp(joint(k).type,'T'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);            
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dij =  r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            ss_omega_i = scew_sym(omega{i});
            ss_omega_j = scew_sym(omega{j});
            % === find normal vectors to hi
            fi = (joint(k).r1 - joint(k).p1)';
            gi = -(scew_sym(fi)*hi);
            % === find normal vector to hj
            fj = (joint(k).r2 - joint(k).p2)';
            % ====================== EQN 1
            gamma(m+1,1) = -hj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*fi + 2*omega{j}'*scew_sym(hj)*A{j}'*A{i}*scew_sym(fi)*omega{i};
            % ====================== EQN 2
            gamma(m+2,1) = -hj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*gi + 2*omega{j}'*scew_sym(hj)*A{j}'*A{i}*scew_sym(gi)*omega{i};
            % ====================== EQN 3
            tmp1 = 2*omega{i}'*scew_sym(fi)*A{i}'*(r_prime{i}-r_prime{j}) + 2*s_jP'*ss_omega_j*A{j}'*A{i}*fi;
            tmp2 = -s_iP'*ss_omega_i*ss_omega_i*fi - s_jP'*ss_omega_j*ss_omega_j*A{j}'*A{i}*fi - dij'*A{i}*ss_omega_i*ss_omega_i*fi;
            gamma(m+3,1) = tmp1 + tmp2;
            clear tmp1 tmp2
            % ====================== EQN 4
            tmp1 = 2*omega{i}'*scew_sym(gi)*A{i}'*(r_prime{i}-r_prime{j}) + 2*s_jP'*ss_omega_j*A{j}'*A{i}*gi;
            tmp2 = -s_iP'*ss_omega_i*ss_omega_i*gi - s_jP'*ss_omega_j*ss_omega_j*A{j}'*A{i}*gi - dij'*A{i}*ss_omega_i*ss_omega_i*gi;
            gamma(m+4,1) = tmp1 + tmp2;
            clear tmp1 tmp2
            % ====================== EQN 5
            gamma(m+5,1) = -fj'*(A{j}'*A{i}*ss_omega_i*ss_omega_i + ss_omega_j*ss_omega_j*A{j}'*A{i})*fi + 2*omega{j}'*scew_sym(fj)*A{j}'*A{i}*scew_sym(fi)*omega{i};
            m=m+5;
            clear i j hi hj s_iP s_jP dij ss_omega_i ss_omega_j fi gi fj
        end
        if(strcmp(joint(k).type,'DI'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);            
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dij =  r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            C = joint(k).dist;
            ss_omega_i = scew_sym(omega{i});
            ss_omega_j = scew_sym(omega{j});
            % ====================== EQN 1
            tmp1 = -2*(r_prime{j}-r_prime{i})'*(r_prime{j}-r_prime{i}) + 2*s_jP'*ss_omega_j*ss_omega_j*s_jP;
            tmp2 = 2*s_iP'*ss_omega_i*ss_omega_i*s_iP - 4*s_jP'*ss_omega_j*A{j}'*A{i}*ss_omega_i*s_iP;
            tmp3 = 4*(r_prime{j}-r_prime{i})'*(A{j}*scew_sym(s_jP)*omega{j} - A{i}*scew_sym(s_iP)*omega{i});
            tmp4 = -2*dij'*(A{i}*ss_omega_i*scew_sym(s_iP)*omega{i} - A{j}*ss_omega_i*scew_sym(s_jP)*omega{j});
            gamma(m+1,1) = tmp1 + tmp2 + tmp3 + tmp4;
            m=m+1;
        end
        if(strcmp(joint(k).type,'G'))
            i=joint(k).bodies(1,1);
            ss_omega_i = scew_sym(omega{i});
            % ====================== EQN 1-2-3
            gamma(m+1:m+3,1) = -A{i}*ss_omega_i*ss_omega_i*joint(k).p1';
            % ====================== EQN 4-5-6
            gamma(m+4:m+6,1) = -0.5*(scew_sym([p_prime{i}(2,1); p_prime{i}(3,1); p_prime{i}(4,1)]) + p_prime{i}(1,1)*eye(3))*omega{i};
            m=m+6;
            clear i ss_omega_i
        end
        if(strcmp(joint(k).type,'RDRV'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            hi = joint(k).q1;
            ss_omega_i = scew_sym(omega{i});
            ss_omega_j = scew_sym(omega{j}); 
            % ====================== EQN 1
            gamma(m+1,1) = -hi*(A{i}'*A{j}*ss_omega_j - ss_omega_i*A{i}'*A{j})*omega{j};            
            m=m+1;
            clear i j hi ss_omega_i ss_omega_j
        end        
        if(strcmp(joint(k).type,'ABX'))
            i=joint(k).bodies(1,1);
            s_iP = [joint(k).dist;0;0];
            ss_omega_i = scew_sym(omega{i});
            tmp =  -A{i}*ss_omega_i*ss_omega_i*joint(k).p1';
            gamma(m+1,1) = tmp(1,1);
            clear i s_iP ss_omega_i tmp
            m=m+1;            
        end
        if(strcmp(joint(k).type,'ABY'))
            i=joint(k).bodies(1,1);
            s_iP = [joint(k).dist;0;0];
            ss_omega_i = scew_sym(omega{i});
            tmp =  -A{i}*ss_omega_i*ss_omega_i*joint(k).p1';
            gamma(m+1,1) = tmp(2,1);
            clear i s_iP ss_omega_i tmp
            m=m+1;          
        end
        if(strcmp(joint(k).type,'ABZ'))
            i=joint(k).bodies(1,1);
            s_iP = [joint(k).dist;0;0];
            ss_omega_i = scew_sym(omega{i});
            tmp =  -A{i}*ss_omega_i*ss_omega_i*joint(k).p1';
            gamma(m+1,1) = tmp(3,1);
            clear i s_iP ss_omega_i tmp
            m=m+1;          
        end
        if(strcmp(joint(k).type,'ANGFIX'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            gamma(m+1:m+3,1) = zeros(3,1); 
            clear i j
            m=m+3;            
        end        
    end
    %% add gamma for math constraints on quaternions
    for i=1:nBody
        gamma(m+i,1) = -2.0*p_prime{i}'*p_prime{i};
    end
end