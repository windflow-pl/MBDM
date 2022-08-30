function Phi_t= calc_Phi_t(joint,nJoint,nCon)
    % allocate
    Phi_t = zeros(nCon,1);
    m=0;
    for k=1:nJoint
        if(strcmp(joint(k).type,'R'))
            m=m+5;
        end
        if(strcmp(joint(k).type,'C'))
            m=m+4;
        end
        if(strcmp(joint(k).type,'S'))
            m=m+3;
        end
        if(strcmp(joint(k).type,'RC'))
            m=m+3;
        end
        if(strcmp(joint(k).type,'T'))
            m=m+5;
        end
        if(strcmp(joint(k).type,'DI'))
            m=m+1;
        end
        if(strcmp(joint(k).type,'G'))
            m=m+6;
        end
        if(strcmp(joint(k).type,'RDRV'))
            Phi_t(m+1,1) = -joint(k).dist;
            m=m+1;
        end
        if(strcmp(joint(k).type,'ABX'))
            m=m+1;
        end
        if(strcmp(joint(k).type,'ABY'))
            m=m+1;
        end
        if(strcmp(joint(k).type,'ABZ'))
            m=m+1;
        end        
    end
end