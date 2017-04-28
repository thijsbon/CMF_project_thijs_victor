%% Velocity profile plot
figure(1)
hold on
plot(u,zc)
grid on

%% semilogx plot
figure(2)
semilogx(y_plus_D(2:end),u(2:end));