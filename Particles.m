rng(1)
clearvars -except l_effective u mu nu Delta_t Steady_State_on zc Time_steps Nz
% Properties of Particles
Np = 100;                 % Number of Particles
Dp = 10e-6;                 % Diameter of Particle
Vp = 4/3*pi*(Dp/2)^3;       % Volume particle
rho_p = 20000;           % Density of particle
C_stokes = 3*pi*mean(mu)*Dp;           % TO BE CHANGED
g = 9.81;               % Gravitational acceleration
Time_steps_2 = Time_steps.*(Steady_State_on==0) + (Steady_State_on==1);
% INSERT (DUST) PARTICLES
Xp = zeros(Np,Time_steps_2+1);
Yp = Xp; 
Zp = Xp;
Xp(:,1) = 0 + randn(Np,1);   % X start position of particles
Yp(:,1) = 0 + randn(Np,1);   % Y start position of particles
Zp(:,1) = 50 + randn(Np,1);   % Z start position of particles
Vpx = zeros(Np,Time_steps_2+1);   % X velocity
Vpy = Vpx;                      % Y velocity
Vpz = Vpx;                      % Z velocity
Vpz(:,1) = 10 + 0.01*randn(Np,1);
a_x = zeros(Np,Time_steps_2+1);   % X acceleration
a_y = a_x;                      % Y
a_z = a_x;                      % Z
F_stokes_x = a_x;
F_stokes_y = a_x;
F_stokes_z = a_x;
F_gravity = a_x;
distance_to_grid = zeros(Np,Time_steps_2);

Time = zeros(Time_steps_2+1,1);     % TIME DUMMY
Euler_Lagrangian_Eddy = 2;          % 1 voor Langevin equation 2 voor discrete eddy model
if Euler_Lagrangian_Eddy == 1 
     %% Langevin model
     V_t = a_x;
     residence_time_x = V_t; residence_time_y = V_t; residence_time_z = V_t; a_eddy_x = V_t; a_eddy_y = V_t; a_eddy_z = V_t;
     u_prime = Vpx; v_prime = Vpx; w_prime = Vpx; b_eddy = a_eddy_x; Ti_star = zeros(Np,Time_steps_2);eddy_life_time = V_t; 
elseif Euler_Lagrangian_Eddy == 2
    V_t = a_x;
    c_t = 0.3;
    T_I_x = zeros(Np,Time_steps_2); T_I_y = T_I_x; T_I_z = T_I_y;
    residence_time_x = V_t; residence_time_y = V_t; residence_time_z = V_t;
    Eddy_start_x = zeros(Np,1);
    Eddy_start_y = zeros(Np,1);
    Eddy_start_z = zeros(Np,1);
    u_prime = Vpx; v_prime = Vpx; w_prime = Vpx;
end
u_prime = 0.01*randn(Np,1);
v_prime = 0.01*randn(Np,1);
w_prime = 0.001*randn(Np,1);

% Time_steps_2 = 2/Delta_t;
for t=1:Time_steps_2
    for particle=1:Np
        %distance_to_grid(particle,t) = find(min(abs(Zp(particle,t)-zc))==abs(Zp(particle,t)-zc)).*(Zp(particle,t)>0) + (Zp(particle,t)<0);
        distance_to_grid(particle,t) = max(find(min(abs(Zp(particle,t)-zc))==abs(Zp(particle,t)-zc)));
        if distance_to_grid(particle,t)==Nz+2
            distance_to_grid(particle,t)=Nz+1;
        end
    end
    %distance_to_grid(:,t) = max(distance_to_grid(:,t));
    if Euler_Lagrangian_Eddy == 1
       %% Langevin model
       V_t(:,t) = (Zp(:,t)<1000).*l_effective(distance_to_grid(:,t),t).*(u(distance_to_grid(:,t)+1,t)-u(distance_to_grid(:,t)-1,t))./zc(distance_to_grid(:,t))' + 1e5*(Zp(:,t)>1000);
       eddy_life_time(:,t)  = c_t * l_effective(distance_to_grid(:,t))./V_t(:,t);
       residence_time_x(:,t)  = l_effective(distance_to_grid(:,t),t)./(abs(Vpx(:,t)-u(distance_to_grid(:,t))));
       residence_time_y(:,t)  = l_effective(distance_to_grid(:,t),t)./(abs(Vpy(:,t)-u(distance_to_grid(:,t))));
       residence_time_z(:,t)  = l_effective(distance_to_grid(:,t))./(abs(Vpz(:,t)-u(distance_to_grid(:,t))));
       Ti_star(:,t) = l_effective(distance_to_grid(:,t))./V_t(:,t);
       a_eddy_x(:,t) = exp(-Delta_t./Ti_star(:,t));
       a_eddy_y(:,t) = exp(-Delta_t./Ti_star(:,t));
       a_eddy_z(:,t) = exp(-Delta_t./Ti_star(:,t));
       %b_eddy(:,t) = 1./(u(distance_to_grid(:,t)+1)-u(distance_to_grid(:,t)-1)./zc(distance_to_grid(:,t))).*sqrt(1-mean([a_eddy_x(:,t);a_eddy_y(:,t);a_eddy_z(:,t)]).^2);
       b_eddy(:,t) = sqrt(2000)*V_t(:,t).*sqrt(1-mean([a_eddy_x(:,t);a_eddy_y(:,t);a_eddy_z(:,t)]).^2);
       u_prime(:,t+1) = a_eddy_x(:,t).*u_prime(:,t)+b_eddy(:,t).*randn(Np,1);
       v_prime(:,t+1) = a_eddy_y(:,t).*v_prime(:,t)+b_eddy(:,t).*randn(Np,1);
       w_prime(:,t+1) = a_eddy_z(:,t).*w_prime(:,t)+b_eddy(:,t).*randn(Np,1);
       %if w_prime(:,t+1)>2*max(max(u))
       %    w_prime(:,t+1) = w_prime(:,t+1)/10;
       %elseif w_prime(:,t+1)<-2*max(max(u))
       %    w_prime(:,t+1) = w_prime(:,t+1)/10;
       %end
    elseif Euler_Lagrangian_Eddy == 2
       V_t(:,t) = l_effective(distance_to_grid(:,t),t).*(u(distance_to_grid(:,t)+1,t)-u(distance_to_grid(:,t)-1,t))./zc(distance_to_grid(:,t))';
       eddy_life_time(:,t)  = c_t * l_effective(distance_to_grid(:,t))./V_t(:,t);
       residence_time_x(:,t)  = l_effective(distance_to_grid(:,t),t)./(abs(Vpx(:,t)-u(distance_to_grid(:,t))));
       residence_time_y(:,t)  = l_effective(distance_to_grid(:,t),t)./(abs(Vpy(:,t)-u(distance_to_grid(:,t))));
       residence_time_z(:,t)  = l_effective(distance_to_grid(:,t))./(abs(Vpz(:,t)-u(distance_to_grid(:,t))));
       T_I_x(:,t) = min(eddy_life_time(:,t),residence_time_x(:,t));
       T_I_y(:,t) = min(eddy_life_time(:,t),residence_time_y(:,t));
       T_I_z(:,t) = min(0.1,min(eddy_life_time(:,t),residence_time_z(:,t)));
       u_prime(:,t+1) = u_prime(:,t).*(t*Delta_t<Eddy_start_x+T_I_x(:,t)) + randn(Np,1).*(t*Delta_t>Eddy_start_x+T_I_x(:,t));
       v_prime(:,t+1) = v_prime(:,t).*(t*Delta_t<Eddy_start_y+T_I_y(:,t)) + randn(Np,1).*(t*Delta_t>Eddy_start_y+T_I_y(:,t));
       w_prime(:,t+1) = w_prime(:,t).*(t*Delta_t<Eddy_start_z+T_I_z(:,t)) + 1.5*randn(Np,1).*(t*Delta_t>Eddy_start_z+T_I_z(:,t));
       % New Eddy start times
       Eddy_start_x = (t*Delta_t>Eddy_start_x+T_I_x(:,t)).*t*Delta_t;
       Eddy_start_y = (t*Delta_t>Eddy_start_y+T_I_x(:,t)).*t*Delta_t;
       Eddy_start_z = (t*Delta_t>Eddy_start_z+T_I_x(:,t)).*t*Delta_t;
    end
    
    
    
    % Calculate forces
    F_stokes_x(:,t) = C_stokes * (u(distance_to_grid(:,t),t)+u_prime(:,t+1)-Vpx(:,t));
    F_stokes_y(:,t) = C_stokes * (v_prime(:,t+1)-Vpy(:,t));
    F_stokes_z(:,t) = C_stokes * (w_prime(:,t+1)-Vpz(:,t));
    F_gravity(:,t) = -g*rho_p*Vp;
    F_stokes_x(:,t) = (Zp(:,t)>0).*F_stokes_x(:,t);
    F_stokes_y(:,t) = (Zp(:,t)>0).*F_stokes_y(:,t);
    F_stokes_z(:,t) =(Zp(:,t)>0).*F_stokes_z(:,t);
    F_gravity(:,t) = (Zp(:,t)>0).*F_gravity(:,t);
    % Calculate accelerations 
    a_x(:,t) = F_stokes_x(:,t)/(Vp*rho_p);
    a_y(:,t) = F_stokes_y(:,t)/(Vp*rho_p);
    a_z(:,t) = (F_stokes_z(:,t)+F_gravity(:,t))/(Vp*rho_p);
    % Calculate velocities
    Vpx(:,t+1) = Vpx(:,t) + a_x(:,t)*Delta_t;
    Vpy(:,t+1) = Vpy(:,t) + a_y(:,t)*Delta_t;
    Vpz(:,t+1) = Vpz(:,t) + a_z(:,t)*Delta_t;
    
    % Calculate new positions
    Xp(:,t+1) = Xp(:,t) + Vpx(:,t+1)*Delta_t + 0.5*a_x(:,t)*Delta_t^2;
    Yp(:,t+1) = Yp(:,t) + Vpy(:,t+1)*Delta_t + 0.5*a_y(:,t)*Delta_t^2;
    Zp(:,t+1) = Zp(:,t) + Vpz(:,t+1)*Delta_t + 0.5*a_z(:,t)*Delta_t^2;
    % Check if particle does not go below ground
    Zp(:,t+1) = (Zp(:,t+1)>0).*Zp(:,t+1);
    Vpz(:,t+1) = (Zp(:,t+1)>0).*Vpz(:,t+1);
    
    % Check if particle does not go above 1000m
    w_prime(:,t) = (Zp(:,t+1)<1000).*w_prime(:,t);
    
    Time(t+1) = t*Delta_t;
end