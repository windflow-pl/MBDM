function [T1,T2,T3,an1,an2,an3] = calc_Thrust2(t)
%     a0 = 5.471261419731895e+02;
%     a1 = -9.224386003387494e-01;
%     b1 = 7.481387825666774e+00;
%     a2 = 5.687687399803831e+00;
%     b2 = 2.838357616988216e-02;
%     a3 = -2.206893100889292e+00;
%     b3 = -3.708819643504395e+00;
%     a4 = -1.725334280068869e+00;
%     b4 = 2.333835354164476e+00;
%     a5 = 1.677366004171450e+00;
%     b5 = 6.277033100479327e-01;
%     w = 2.445983032762768e-02;
%     % === rotational speed [rpm]
%     rpm = 12.1;
%     omega = rpm*2*pi/60; %[rad/s]
%     % === blade 1 [0 angle]
%     x0=0.0;%[rad]
%     x=(x0+omega*t)*360/2/pi;%[deg]
%     while (x>360)
%         x = x-360;
%     end
%     an1 = x;
%     T1  =  156.934283047*(a0+a1*cos(x*w)+b1*sin(x*w)+a2*cos(2*x*w)+b2*sin(2*x*w)+a3*cos(3*x*w)+b3*sin(3*x*w)+a4*cos(4*x*w)+b4*sin(4*x*w)+a5*cos(5*x*w)+b5*sin(5*x*w));
%     % === blade 2 [240 angle]
%     x0=240/360*2*pi;%[rad]
%     x=(x0+omega*t)*360/2/pi;%[deg]
%     while (x>360)
%         x = x-360;
%     end
%     an2 = x;
%     T2  =  156.934283047*(a0+a1*cos(x*w)+b1*sin(x*w)+a2*cos(2*x*w)+b2*sin(2*x*w)+a3*cos(3*x*w)+b3*sin(3*x*w)+a4*cos(4*x*w)+b4*sin(4*x*w)+a5*cos(5*x*w)+b5*sin(5*x*w));
%     % === blade 3 [120 angle]
%     x0=120/360*2*pi;%[rad]
%     x=(x0+omega*t)*360/2/pi;%[deg]
%     while (x>360)
%         x = x-360;
%     end
%     an3 = x;
%     T3  =  156.934283047*(a0+a1*cos(x*w)+b1*sin(x*w)+a2*cos(2*x*w)+b2*sin(2*x*w)+a3*cos(3*x*w)+b3*sin(3*x*w)+a4*cos(4*x*w)+b4*sin(4*x*w)+a5*cos(5*x*w)+b5*sin(5*x*w));
    % ======== THRUST
%     Thrust = [T1+T2+T3;0;0];

%     T1 = 260000;
%     T2 = 260000;
%     T3 = 260000;
%     an1 = 0;
%     an2 = 120;
%     an3 = 240;

    % === average value
    T1 = 8.604254245799026e+04; 
    T2 = 8.604254245799026e+04; 
    T3 = 8.604254245799026e+04; 
    an1 = 0;
    an2 = 120;
    an3 = 240;
end



