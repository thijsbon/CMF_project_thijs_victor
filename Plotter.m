%% Velocity profile plot
figure(1)
hold on

plot(u(2:end-1,end),zc(2:end-1))
grid on

%% semilogx plot
figure(2)
semilogx(y_plus_D(2:end),u(2:end,end));

plot(u(2:end-1,end),zc(2:end-1),'ob');
xlabel('u(m/s)'); ylabel('z(m)');
grid on
%ANALYTIC SOLUTION
if show_analytic==1
    u_analytic = -6*Q/H^3*(zc.^2-H*zc)*(prescribeswitch == 1)*(bcswitch == 0)*(uwall1 == 0) + ...
             1./(2.*mu).*dpdx.*(zc.^2-H*zc)*(prescribeswitch == 0)*(bcswitch == 0)*(uwall1 == 0); %poiseuille flow, i.e. v = 0 at walls, if Q is specified

    plot(u_analytic(2:end-1),zc(2:end-1),'r'); legend('numerical','analytic');
    figure(2)
    hold on
    plot(u(2:end-1)/max(u),zc(2:end-1),'ob')
    plot(u_analytic(2:end-1)/max(u_analytic),zc(2:end-1),'r')
    legend('numerical','analytic')
    xlabel('u/u_{max}'); ylabel('z(m)')
    grid on
    figure(3)
    hold on
    plot(u(2:end-1,end)/mean(mean(u)),zc(2:end-1),'ob')
    plot(u_analytic(2:end-1)/mean(u_analytic),zc(2:end-1),'r')
    legend('numerical','analytic')
    xlabel('u/<u>'); ylabel('z(m)')
    grid on
end

% %%
figure(4)
semilogx(zc(2:end-1),u(2:end-1,end),'ob');
%%
figure(5)
plot(zf(2:end-1),nu_t,'o')

%% Transient plotter
if Steady_State_on == 0 
    amount_of_lines = 19;
    colour_1 = 0:1/Time_steps:1;
    colour_2 = [colour_1' flip(colour_1') flip(abs(-colour_1'))];
    figure(6)
    hold on
    i=1;
    for col = 1:round(Time_steps/(amount_of_lines-1)):Time_steps        
        plot(u(2:end-1,col),zc(2:end-1),'Color',colour_2(col,:))
        legendInfo{i} = ['t = ' num2str((col-1)*Delta_t) 's'];
        %legend_string(col) = t = , (col+1)*Delta_t);
        i = i+1;
    end
    legend(legendInfo)
    title(['Changing Velocity profile for \omega = ',num2str(omega_unsteady) ' Rad/s'])
    xlabel('Velocity (m/s)')
    ylabel('Channel Height (m)')
end