%% Velocity profile plot
figure(1)
hold on
plot(u(2:end-1,end),zc(2:end-1))
grid on

%% semilogx plot
figure(2)
semilogx(y_plus_D(2:end),u(2:end,end));