rng(1)
% Properties of Particles
Np = 1;                 % Number of Particles
Dp = 1e-5;                 % Diameter of Particle
Vp = 4/3*pi*(Dp/2)^3;       % Volume particle
rho_p = 20000;           % Density of particle
C_stokes = 3*pi*mean(mu)*Dp;           % TO BE CHANGED
g = 9.81;               % Gravitational acceleration
% INSERT (DUST) PARTICLES
Xp = zeros(Np,Time_steps+1);
Yp = Xp; 
Zp = Xp;
Xp(:,1) = 0 + randn(Np,1);   % X start position of particles
Yp(:,1) = 0 + randn(Np,1);   % Y start position of particles
Zp(:,1) = 20 + randn(Np,1);   % Z start position of particles
Vpx = zeros(Np,Time_steps+1);   % X velocity
Vpy = Vpx;                      % Y velocity
Vpz = Vpx;                      % Z velocity
a_x = zeros(Np,Time_steps+1);   % X acceleration
a_y = a_x;                      % Y
a_z = a_x;                      % Z

Time = zeros(Time_steps+1,1);     % TIME DUMMY
for t=1:Time_steps
    distance_to_grid = find(min(abs(Zp(:,t)-zc))==abs(Zp(:,t)-zc));
    % Calculate forces
    F_stokes_x(:,t) = C_stokes * (u(distance_to_grid(end),t)-Vpx(:,t));
    F_stokes_y(:,t) = C_stokes * -Vpy(:,t);
    F_stokes_z(:,t) = C_stokes * -Vpz(:,t);
    F_gravity(:,t) = -g*rho_p*Vp;
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
    
    Time(t+1) = t*Delta_t;
end