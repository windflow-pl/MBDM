function A = calc_A_angle(angle)
% function calculate rotational matrix based on 1-2-3
% euler angles
    Ax = [ 1           , 0             ,              0;
           0           , cos(angle(1)) , -sin(angle(1));                 
           0           , sin(angle(1)) ,  cos(angle(1))];               
    Ay = [cos(angle(2)), 0             , -sin(angle(2));
          0            , 1             ,              0;
          sin(angle(2)), 0             ,  cos(angle(2))];
    Az = [cos(angle(3)), -sin(angle(3)),              0;
          sin(angle(3)), cos(angle(3)) ,              0;
          0            , 0             ,              1];
    A = Ax*Ay*Az;  
end