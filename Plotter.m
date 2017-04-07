%% Velocity profile plot
figure(1)
hold on
plot(u(2:end-1),zc(2:end-1),'ob');
xlabel('u(m/s)'); ylabel('z(m)');
grid on
%ANALYTIC SOLUTION
u_analytic = -6*Q/H^3*(zc.^2-H*zc)*(prescribeswitch == 1)*(bcswitch == 0)*(uwall1 == 0) + ...
             1./(2.*mu).*dpdx.*(zc.^2-H*zc)*(prescribeswitch == 0)*(bcswitch == 0)*(uwall1 == 0); %poiseuille flow, i.e. v = 0 at walls, if Q is specified
if show_analytic==1
plot(u_analytic(2:end-1),zc(2:end-1),'r'); legend('numerical','analytic');

end
figure(2)
hold on
plot(u(2:end-1)/max(u),zc(2:end-1),'ob')
plot(u_analytic(2:end-1)/max(u_analytic),zc(2:end-1),'r')
legend('numerical','analytic')
xlabel('u/u_{max}'); ylabel('z(m)')
grid on

% %%
figure
semilogx(zc(2:end-1),u(2:end-1),'ob');
%%
figure
hold on
plot(u(2:end-1)/mean(u),zc(2:end-1),'ob')
plot(u_analytic(2:end-1)/mean(u_analytic),zc(2:end-1),'r')
legend('numerical','analytic')
xlabel('u/<u>'); ylabel('z(m)')
grid on
%%
figure
plot(zf(2:end-1),nu_t,'o')
