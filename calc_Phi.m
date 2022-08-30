function [Phi, er] = calc_Phi(joint,q,t,nBody,nJoint,nCon)
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
    % alloc A matrices
    for i=1:nBody
        A{i} = zeros(3,3); 
        A{i} = calc_A(p{i});
    end
    % allocate Phi
    Phi = zeros(nCon,1);
    % set error handler to 0 => no error
    er = 0;
    % calc Phi
    m =0;
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
            if(gi'*hi~=0 || fi'*hi~=0 )
                er = 1;
                fprintf('Did not properly find normal vectors to vector hi in calc_Phi for joint %d !!!\n',k);
                return;
            end
            % === UPDATE PHI
            Phi(m+1:m+3,1) = r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            Phi(m+4,1) = fi'*A{i}'*A{j}*hj;
            Phi(m+5,1) = gi'*A{i}'*A{j}*hj;
            m=m+5;
            clear i j hi hj s_jP  s_iP fi gi
        end 
        if(strcmp(joint(k).type,'C'))            
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
            if(gi'*hi~=0 || fi'*hi~=0 )
                er = 1;
                fprintf('Did not properly find normal vectors to vector hi in calc_Phi for joint %d !!!\n',k);
                return;
            end
            % === UPDATE PHI
            Phi(m+1,1) = fi'*A{i}'*A{j}*hj;
            Phi(m+2,1) = gi'*A{i}'*A{j}*hj;
            Phi(m+3,1) = fi'*A{i}'*dij;
            Phi(m+4,1) = gi'*A{i}'*dij;
            m=m+4;
            clear i j hi hj s_jP  s_iP fi gi
        end 
        if(strcmp(joint(k).type,'S'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            % === UPDATE PHI
            Phi(m+1:m+3,1) = r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            m=m+3;
            clear i j s_iP s_jP
        end
        if(strcmp(joint(k).type,'RC'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            hi = (joint(k).q1 - joint(k).p1)';
            hj = (joint(k).q2 - joint(k).p2)';
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dji =  r{i} + A{i}*s_iP - r{j} - A{j}*s_jP;
            % === find normal vectors to hj
            fj = (joint(k).r2 - joint(k).p2)';
            gj = -(scew_sym(fj)*hj);
            % === UPDATE PHI
            Phi(m+1,1) = hi'*A{i}'*A{j}*hj;
            Phi(m+2,1) = fj'*A{j}'*dji;
            Phi(m+3,1) = gj'*A{j}'*dji;
            clear i j hi hj s_jP s_iP fj gj dji
            m=m+3;
        end
        if(strcmp(joint(k).type,'T'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
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
            % === UPDATE PHI
            Phi(m+1,1) = fi'*A{i}'*A{j}*hj;            
            Phi(m+2,1) = gi'*A{i}'*A{j}*hj;            
            Phi(m+3,1) = fi'*A{i}'*dij;
            Phi(m+4,1) = gi'*A{i}'*dij;
            Phi(m+5,1) =  fi'*A{i}'*A{j}*fj;
            clear i j hi hj s_jP s_iP fj fi gi dij
            m=m+5;
        end
        if(strcmp(joint(k).type,'DI'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            s_iP = joint(k).p1';
            s_jP = joint(k).p2';
            dij =  r{j} + A{j}*s_jP - r{i} - A{i}*s_iP;
            C = joint(k).dist;
            % === UPDATE PHI
            Phi(m+1,1) = dij'*dij - C*C;
            m=m+1;
            clear i j s_jP s_iP dij C
        end
        if(strcmp(joint(k).type,'G'))
            i=joint(k).bodies(1,1);
            clear j
            Phi(m+1:m+3,1) = r{i} - joint(k).p1';
            Phi(m+4:m+6,1) = [p{i}(2,1);p{i}(3,1);p{i}(4,1)] - joint(k).q1';
            m=m+6;
            clear i
        end
        if(strcmp(joint(k).type,'RDRV'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            fj = joint(k).p2';
            fi = joint(k).p1';
            hi = joint(k).q1';
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
            % == find angle
            if (s>=0 && c>=0)
                ang = asin(s);
            elseif (s>=0 && c<0)
                ang = pi - asin(s);
            elseif (s<0 && c<0)
                ang = pi - asin(s);
            elseif (s<0 && c>=0)
                ang = 2*pi + asin(s);
            end
            % check angle to be real number            
            if (~isreal(ang))
                er = 1;
            end
            % == account for full revolutions
            n=0;
            C_fun = joint(k).r1(1,1) + joint(k).dist*t; 
            % == if rotating counter clockwise
            for tmp=1:100
                if (C_fun-2*n*pi>=0 && C_fun-2*n*pi<2*pi)
                    break;
                else
                    % == if rotating counter clockwise
                    if (joint(k).dist>0)
                        n=n+1;
                    end
                    % == if rotating clockwise
                    if (joint(k).dist<0)
                        n=n-1;
                    end
                end
            end 
            % == help N-R to converge
            % calculate Phi
            if ((ang - C_fun + 2*n*pi)==2*pi)
                Phi(m+1,1) = 0.0;
            else
                Phi(m+1,1) = ang - C_fun + 2*n*pi;
            end 
            % a bit more help, by choosing smaller angle to correct
            if(Phi(m+1,1)<-pi)
                Phi(m+1,1) = 2*pi + (ang - (C_fun) + 2*n*pi);
            end
            if(Phi(m+1,1)>pi)
                Phi(m+1,1) = 2*pi - (ang - (C_fun) + 2*n*pi);
            end 
            m=m+1;
            clear i j fj fi gi hi c s ang C_fun n
        end
        if(strcmp(joint(k).type,'ABX'))
            i=joint(k).bodies(1,1);
            Phi(m+1,1) = r{i}(1,1) - joint(k).dist;
            m=m+1;
            clear i
        end
        if(strcmp(joint(k).type,'ABY'))
            i=joint(k).bodies(1,1);
            Phi(m+1,1) = r{i}(2,1) - joint(k).dist;
            m=m+1;
            clear i
        end
        if(strcmp(joint(k).type,'ABZ'))
            i=joint(k).bodies(1,1);
            Phi(m+1,1) = r{i}(3,1) - joint(k).dist;
            m=m+1;
            clear i
        end
        if(strcmp(joint(k).type,'ANGFIX'))
            i=joint(k).bodies(1,1);
            j=joint(k).bodies(1,2);
            % 1 denotes orientation of master body
            % 2 denotes additional rotation relative to the master
            e01 = joint(k).q2(1);
            e11 = joint(k).p2(1);
            e21 = joint(k).p2(2);
            e31 = joint(k).p2(3);
            e02 = joint(k).q2(2);
            e12 = joint(k).r2(1);
            e22 = joint(k).r2(2);
            e32 = joint(k).r2(3);
                        
            % elements of resulting quaternion
            el1 = e22*e31 - e32*e21 + e02*e11 + e01*e12;
            el2 = e32*e11 - e12*e31 + e02*e21 + e01*e22;
            el3 = e12*e21 - e22*e11 + e02*e31 + e01*e32;
            
            Phi(m+1:m+3,1) = joint(k).p1'-[el1; el2; el3];
            m=m+3;
            clear i j
        end
        
    end
    %% add math constraints
    for i=1:nBody
        Phi(m+1,1) = p{i}'*p{i} -1;
        m = m+1;
    end

end