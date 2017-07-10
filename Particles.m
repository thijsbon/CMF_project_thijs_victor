rng(1)
clearvars -except l_effective u mu nu Delta_t Steady_State_on zc Time_steps Nz rho dz
%%
% Properties of Particles
Np = 10;                 % Number of Particles
np = 1000;                % 1 particle represents np particles
Dp = 1e-4;                 % Diameter of Particle
Vp = 4/3*pi*(Dp/2)^3;       % Volume particle
rho_p = 2000;           % Density of particle
mu = mean(mu);
%%
rain_on = 1;        % 1 on, 0 off
Nd = 10;                %Number of druppels
nd = 1000;              %1Nd Represents nd particles
Dd = 1e-3;              %Diameter of druppel, range 5e-4 -> 5e-3
Vd = 4/3*pi*(Dd/2)^3;   %Volume of druppel
rho_d = 1000;           %Density of druppel
C_stokes_rain = 3*pi*mu*Dd; %Stokes drag coefficient for droplets

%C_stokes = 3*pi*mu*Dp^2;           % TO BE CHANGED
C_stokes = 3*pi*mu*Dp;
g = 9.81;               % Gravitational acceleration
Time_steps_2 = Time_steps.*(Steady_State_on==0) + (Steady_State_on==1);
Delta_Time_for_particles = 0.001;
Time_steps_for_particles = round(Delta_t/Delta_Time_for_particles);
% Particle-Particle collisions
Dp_effective = sqrt(np)*Dp;
Dd_effective = sqrt(nd)*Dd;

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
Vpz(:,1) = 0 + 0.01*randn(Np,1);
Vpx(:,1) = u(77,1)+randn(Np,1);
a_x = zeros(Np,Time_steps_2+1);   % X acceleration
a_y = a_x;                      % Y
a_z = a_x;                      % Z
F_stokes_x = a_x;
F_stokes_y = a_x;
F_stokes_z = a_x;
F_gravity = a_x;
distance_to_grid = zeros(Np,Time_steps_2);

% INSERT DROPLETS
Xd = zeros(Nd,Time_steps_2+1);
Yd = Xd;Zd = Xd;
Xd(:,1) = rand(Nd,1);   % X start position of particles
Yd(:,1) = rand(Nd,1);   % Y start position of particles
Zd(:,1) = rand(Nd,1);   % Z start position of particles
Vdx = zeros(Nd,Time_steps_2+1);
Vdy = Vdx;
Vdz = Vdx;
Vdz(:,1) = -2/9*1000/1.5e-5*g*(Dd/2)^2;
ad_x = Vdx;
ad_y = ad_x;
ad_z = ad_x;
F_stokes_dx = ad_x;
F_stokes_dy = ad_x;
F_stokes_dz = ad_x;

Time = zeros(Time_steps_2+1,1);     % TIME DUMMY
Euler_Lagrangian_Eddy = 2;          % 1 voor Langevin equation 2 voor discrete eddy model
u_prime = Vpx; v_prime = Vpx; w_prime = Vpx;
if Euler_Lagrangian_Eddy == 1 
     %% Langevin model
     V_t = a_x;
     c_t = 0.3;
     residence_time_x = V_t; residence_time_y = V_t; residence_time_z = V_t; 
     Ti_star = zeros(Np,Time_steps_2);eddy_life_time = V_t; 
elseif Euler_Lagrangian_Eddy == 2
    V_t = a_x;
    c_t = 0.3;
    eddy_life_time = V_t; 
    T_I_x = zeros(Np,Time_steps_2); T_I_y = T_I_x; T_I_z = T_I_y;
    residence_time_x = V_t; residence_time_y = V_t; residence_time_z = V_t;
    Eddy_start_x = zeros(Np,1);
    Eddy_start_y = zeros(Np,1);
    Eddy_start_z = zeros(Np,1);
end
u_prime(:,1) = randn(Np,1);
v_prime(:,1) = 0.1*randn(Np,1);
w_prime(:,1) = 0.1*randn(Np,1);

%%
upper = 0;
%Time_steps_2 = 10;
for t=1:Time_steps_2
    t
    %% Sub time program
    Delta_Particle;
    %% Update values
    F_stokes_x(:,t+1) = F_stokes_x_tt(:,end);
    F_stokes_y(:,t+1) = F_stokes_x_tt(:,end);
    F_stokes_z(:,t+1) = F_stokes_x_tt(:,end);
    %F_gravity(:,t+1) = F_gravity_tt(:,end);
    a_x(:,t+1) = a_x_tt(:,end);
    a_y(:,t+1) = a_y_tt(:,end);
    a_z(:,t+1) = a_z_tt(:,end);
    Vpx(:,t+1) = Vpx_tt(:,end);
    Vpy(:,t+1) = Vpy_tt(:,end);
    Vpz(:,t+1) = Vpz_tt(:,end);
    Xp(:,t+1) = Xp_tt(:,end);
    Yp(:,t+1) = Yp_tt(:,end);
    Zp(:,t+1) = Zp_tt(:,end);
    u_prime(:,t+1) = u_prime_tt(:,end);
    v_prime(:,t+1) = v_prime_tt(:,end);
    w_prime(:,t+1) = w_prime_tt(:,end);
    eddy_life_time(:,t+1) = eddy_life_time_tt(:,end);
    % Check if particle does not go above 1000m
    %w_prime(:,t) = (Zp(:,t+1)<1000).*w_prime(:,t);
    
     %% Inter particle collisions NIET AF!
%      for particle=1:Np
%          inter_particle_distance(particle,:) = sqrt((Xp(:,t+1)-Xp(particle,t+1)).^2+(Yp(:,t+1)-Yp(particle,t+1)).^2+(Zp(:,t+1)-Zp(particle,t+1)).^2)';
%          %if inter_particle_distance(particle,:)'<100*sqrt((Vpx(:,t+1)-Vpx(particle,t+1)).^2+(Vpy(:,t+1)-Vpy(particle,t+1)).^2+(Vpz(:,t+1)-Vpz(particle,t+1)).^2)*Delta_t
%             probability_ij(particle,:) = np*pi*Dp_effective^2.*sqrt((Vpx(:,t+1)-Vpx(particle,t+1)).^2+(Vpy(:,t+1)-Vpy(particle,t+1)).^2+(Vpz(:,t+1)-Vpz(particle,t+1)).^2)*Delta_t;
%          %end
%      end    
    
    Time(t+1) = t*Delta_t;
end