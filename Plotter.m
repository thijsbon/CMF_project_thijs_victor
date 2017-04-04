%% Velocity profile plot
figure(1)
hold on
plot(u(2:end-1),zc(2:end-1),'ob')
grid on
%ANALYTIC SOLUTION
u_analytic = -6*Q/H^3*(zc.^2-H*zc)*(prescribeswitch == 1)*(bcswitch == 0)*(uwall1 == 0) + ...
             1./(2.*mu).*dpdx.*(zc.^2-H*zc)*(prescribeswitch == 0)*(bcswitch == 0)*(uwall1 == 0); %poiseuille flow, i.e. v = 0 at walls, if Q is specified
if show_analytic==1
plot(u_analytic(2:end-1),zc(2:end-1),'r')
end
figure(2)
hold on
plot(u(2:end-1)/max(u),zc(2:end-1),'ob')
plot(u_analytic(2:end-1)/max(u_analytic),zc(2:end-1),'r')
grid on