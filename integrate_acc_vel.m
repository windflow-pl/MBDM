clc
clear all

src = 'results/data/no_thrust_2.dat';
delimiterIn = ' ';
results = importdata(src,delimiterIn);
% t,x,y,z,vx,vy,vz,ax,ay,az,angx,angy,angz,omegax,omegay,omegaz,epsx,epsy,epsz,V_submerged

a_res = [results(:,8),results(:,9),results(:,10)];
v_res = [results(:,5),results(:,6),results(:,7)];
x_res = [results(:,2),results(:,3),results(:,4)];
dt = results(2,1)-results(1,1);
for i=1:length(a_res);
    v(i,1) = results(i,5) + a_res(i,1)*dt;
    v(i,2) = results(i,6) + a_res(i,2)*dt;
    v(i,3) = results(i,7) + a_res(i,3)*dt;
    x(i,1) = results(i,2) + v_res(i,1)*dt;
    x(i,2) = results(i,3) + v_res(i,2)*dt;
    x(i,3) = results(i,4) + v_res(i,3)*dt;
end

% v = 

figure
hold on
subplot(2,3,1)
plot(results(:,1),v(:,1),'o-',results(:,1),v_res(:,1),'+-')
legend('v_x integrated','v_x MBDM','Location','NorthOutside')

subplot(2,3,2)
plot(results(:,1),v(:,2),'o-',results(:,1),v_res(:,2),'+-')
legend('v_y integrated','v_y MBDM','Location','NorthOutside')

subplot(2,3,3)
plot(results(:,1),v(:,3),'o-',results(:,1),v_res(:,3),'+-')
legend('v_z integrated','v_z MBDM','Location','NorthOutside')

subplot(2,3,4)
plot(results(:,1),x(:,1),'o-',results(:,1),x_res(:,1),'+-')
legend('x integrated','x MBDM','Location','NorthOutside')

subplot(2,3,5)
plot(results(:,1),x(:,2),'o-',results(:,1),x_res(:,2),'+-')
legend('y integrated','y MBDM','Location','NorthOutside')

subplot(2,3,6)
plot(results(:,1),x(:,3),'o-',results(:,1),x_res(:,3),'+-')
legend('z integrated','z MBDM','Location','NorthOutside')


hold off



