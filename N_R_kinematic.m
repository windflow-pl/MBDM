function [q_out,Phi_q,er] = N_R_kinematic(q, joint,nBody,nJoint, nCon, t)
    tolq=1.0e-3;
    tolf=1.0e-3;
    steps = 300;
    relax.a = [0.5,0.25,0.125,0.0625,0.03125];
    relax.step = [50,100,150,200,250];
    relax.index = 1;

    % ALLOCATE
    DPhi = zeros(nCon,7*nBody+1);
    % redirect input to output before going for iteration
    er = 0;
    Phi_q = 0;
    q_out = q;
    for k=1:steps
%         k
        % tmp
        DPhi_1 = DPhi;
        % calc Phi
        [Phi, er] = calc_Phi(joint,q_out,t,nBody,nJoint,nCon);

        %%%%%%%%%%%% tmp
%         Phi
        %%%%%%%%%%%%
        if er==1
           disp(sprintf('ERROR: angle has imaginary part in N-R\n'));
           return;
        end
        % check if Phi[i]<tolf
        j=0;
        for i=1:nCon
            if (abs(Phi(i,1))>=tolf) j=1; end
        end
        if (j==0)
%             disp(sprintf('N-R reached conv in tolf\n'))
            break;
        end
        % place [-Phi] as a last column of DPhi
        for i=1:nCon
            DPhi(i,nCon+1) = -1.0*Phi(i,1);
        end
        % calc Phi_q
        Phi_q = calc_Phi_q(joint,q_out,nBody,nJoint,nCon);
        % place Phi_q in DPhi
        for i=1:nCon
            for j=1:7*nBody
                DPhi(i,j) = Phi_q(i,j);
            end
        end
        % ================ solve
        [A, loc_indx, er] = RowReducedEchelonForm(DPhi,1);
        if er==1
           return;
        end
        dsol = get_solution(A, nCon,7*nBody+1,loc_indx);
%         dsol
        % check if dsol<tolq
        j=0;
        for i=1:nCon
            if(abs(dsol(i)) >=tolq) j=1; end
        end
        if (~isreal(dsol))
            er = 1;
            fprintf('ERROR: dsol has imaginary part in N-R\n');
            return;
        end
        if (j==0)
%             disp(sprintf('N-R reached conv in tolq\n'))
            break;
        end
        % UPDATE
        if k>relax.step(relax.index)
##          k
##          relax.index
##          relax.a(relax.index)
          for i=1:nCon
            q_out(i,1) = q_out(i,1)+dsol(i,1)*relax.a(relax.index);
          end
        else
          for i=1:nCon
              q_out(i,1) = q_out(i,1)+dsol(i,1);
          end
        endif
        % UPDATE relaxation index
        if relax.index<=length(relax.step)
          if k==relax.step(relax.index);
            fprintf('Switching relaxation index to a=%f\n',relax.a(relax.index));
          endif
          if k==relax.step(relax.index+1);
            relax.index = relax.index+1;
            fprintf('Switching relaxation index to a=%f\n',relax.a(relax.index));
          endif
        endif


        if k==steps
            er = 1;
            disp(sprintf('N-R didnt converge in %d steps. Possible solutions: increase ''steps'', decrease time step or check Jakobian.\n',k));
            return;
        end
    end
end
