rng(1)
tic
clearvars -except l_effective u mu nu Delta_t Steady_State_on zc Time_steps Nz rho dz
%% Properties of Particles
Np = 10;                 % Number of Particles
np = 1000;              % 1 particle represents np particles
Dp = 5e-5;              % Diameter of Particle
Vp = 4/3*pi*(Dp/2)^3;   % Volume particle
rho_p = 2000;           % Density of particle
mu = mean(mu);
%% rain
rain_on = 1;        % 1 on, 0 off
Nd = 1;                %Number of druppels
nd = 1000;              %1Nd Represents nd particles
Dd = 1e-3;              %Diameter of druppel, range 5e-4 -> 5e-3
Vd = 4/3*pi*(Dd/2)^3;   %Volume of druppel
rho_d = 1000;           %Density of druppel
C_stokes_rain = 3*pi*mu*Dd; %Stokes drag coefficient for droplets

C_stokes = 3*pi*mu*Dp;
g = 9.81;               % Gravitational acceleration
Time_steps_2 = Time_steps.*(Steady_State_on==0) + (Steady_State_on==1);
Delta_Time_for_particles = 0.001;
if Dp<1e-5
    Delta_Time_for_particles = 0.0002;
end
Time_steps_for_particles = round(Delta_t/Delta_Time_for_particles);

% Speed up, on/off
speed_up = 1;       % 1 = on, 0 = off
                    % When using 1, the sub time is calculated only for 100
                    % time steps and from that a mean velocity and
                    % deviation is calculated, which uses monte carlo to
                    % determine the velocity and location at the next full
                    % time step. USE ONLY WHEN YOU HAVE LIMITED TIME!
                    % Using 0 uses the normal method, disered when you have
                    % the time and resources.
speed_up_method = 7;    % Use 1 for only statistical properties for 
                        % location and velocity. This might have errors
                        % Use 2 to calculate end velocities and locations
                        % with accelerations calculated by statistical
                        % properties. This might be better, however, this
                        % can suffer from extreme accelerations.
                        % Use 3 to use a linear extrapolate model. If
                        % x_end-x_begin < 1*std(x) then the slope is
                        % calculated by a = (x(:,end) - x(:,1)) /
                        % (100*Delta_Time_for_particles -
                        % 1*Delta_Time_for_particles)
                        % And then the final locations are determined using
                        % this extrapolation
                        % If the end values are above/below 1*std, then the
                        % same as setting 1 is used.
                        % the 
                        % Use 4 to do the average of 3 and 1.
                        % Use 5 for terminal velocity model. This model
                        % uses statistics everywhere, except for Zp which
                        % uses zp_tt(:,end)-v_terminal*Delta_t, with
                        % v_terminal = 2/9*Drho/mu *g*(Dp/2)^2
                        % Use 6 for model of 5, where the Zp_tt = Zp_tt +
                        % std(z)
                        % Model 7 is similar to model 5, however, it allows
                        % to increase the average velocity effect with the
                        % add_precentage variable
% LEES DIT!
% Het blijkt uit mijn experimenten dat speed_up_method 1 het goed doet maar
% onderpresteerd. 2 en 3 trekken de particles VEEL te snel naar de grond, 4
% doet het medium, niet ideaal. Nummer 5 doet het erg goed, maar mist de
% fluctuations die je krijgt als je het normaal uit rekent. Daarom doet
% method 6 het heel goed, want die ziet er erg uit zoals het normaal
% uitrekenen, zonder de reken kosten.
% HET IS DAAROM AANGERADEN om model 5 te kiezen, en dan model 6 te runnen.
% Als het verschil tussen 5 en 6 heel groot is, kies nummer 5.
% Als het verschil tussen 5 en 6 klein is, kies nummer 6.
% Je kan altijd spelen met Time_steps_for_particles om het proces te
% versnellen, maar wees bewust van het feit dat de statistics dan fout
% kunnen gaan.
   
add_percentage = 1.06;  % Use 1.06 for a Dp=1e-5
sub_method_2 = 3;       % Only need to use if speed_up_method = 2, if this 
                        % value == 1, then the end location is calculated
                        % according to the right value of Vpx
                        % If the value == 2, the the end location is
                        % calculated with the the last value of Vpx (which
                        % is actually incorrect).
                        % If the value == 3, uses a one time time step

if speed_up == 1
    Original_Time_steps_for_particles = Time_steps_for_particles;
    Time_steps_for_particles = 1000;
    mod_6 = -log(5*Dp); % Constante nodig voor speed_up_method == 6
end

% Particle-Particle collisions
Dp_effective = sqrt(np)*Dp;
Dd_effective = sqrt(nd)*Dd;
probability_ij = zeros(Np,Nd);
V_relative = zeros(Np,Nd);
droplet_particle_distance = zeros(Np,Nd);
dpd_min = zeros(Np,1);
minimum_distance = zeros(Np,1);
mu_d = 0.2;                     % friction coefficient (0 (frictionless sliding) to 0.4 (maximum friction))
e_r = 0.9;                      % restitution coefficient (0.8 to 1 (e_r =1 == frictionless sliding))
Collision = zeros(Np,Time_steps_2);
Chance_to_get_stuck = 10/np;     % 
dtt = Delta_Time_for_particles; % Used delta t in determination of probabilities, influences the amount of collisions very much.

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
Landing_time = ones(Np,2);

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
F_stokes_xd = ad_x;
F_stokes_yd = ad_x;
F_stokes_zd = ad_x;

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
toc
for t=1:Time_steps_2
    tic
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
       
    %% Rain
    if rain_on == 1
        ad_x(:,t+1) = ad_x_tt(:,end);
        ad_y(:,t+1) = ad_x_tt(:,end);
        ad_z(:,t+1) = ad_x_tt(:,end);
        Vdx(:,t+1) = Vdx_tt(:,end);
        Vdy(:,t+1) = Vdy_tt(:,end);
        Vdz(:,t+1) = Vdz_tt(:,end);
        Xd(:,t+1) = Xd_tt(:,end);
        Yd(:,t+1) = Yd_tt(:,end);
        Zd(:,t+1) = Zd_tt(:,end);
        F_stokes_xd(:,t+1) = F_stokes_xd_tt(:,end);
        F_stokes_yd(:,t+1) = F_stokes_yd_tt(:,end);
        F_stokes_zd(:,t+1) = F_stokes_zd_tt(:,end);
    
    
        %% Droplet-Particle Collision
        % Determination of collision
        probability_ij = zeros(Np,Nd);
        for particle=1:Np
            %% Landing time registration (done in this loop to reduce time
            if Zp(particle,t+1)==0 && Zp(particle,t)>0
                Landing_time(particle,1) = (t+1)*Delta_t;
                Landing_time(particle,2) = t+1;
            end
            %% Check if particle doesn't move when it touched the ground
            Xp(particle,t+1) = Xp(particle,t+1).*((Landing_time(particle,2))==1)+Xp(particle,Landing_time(particle,2)).*((Landing_time(particle,2))>1);
            Yp(particle,t+1) = Yp(particle,t+1).*((Landing_time(particle,2))==1)+Yp(particle,Landing_time(particle,2)).*((Landing_time(particle,2))>1);
            Zp(particle,t+1) = Zp(particle,t+1).*((Landing_time(particle,2))==1)+Zp(particle,Landing_time(particle,2)).*((Landing_time(particle,2))>1);

            %% Collision part
            %droplet_particle_distance(particle,:) = sqrt((Xp(particle,t+1)-floor(Xp(particle,t+1))-Xd(:,t+1)).^2+(Yp(particle,t+1)-floor(Yp(particle,t+1))-Yd(:,t+1)).^2+(Zp(particle,t+1)-floor(Zp(particle,t+1))-Zd(:,t+1)).^2);
            droplet_particle_distance(particle,:) = sqrt((Xd(:,t+1)-Xp(particle,t+1)+floor(Xp(particle,t+1))).^2+(Yd(:,t+1)-Yp(particle,t+1)+floor(Yp(particle,t+1))).^2+(Zd(:,t+1)-Zp(particle,t+1)+floor(Zp(particle,t+1))).^2);
            V_relative(particle,:) = (sqrt((Vpx(particle,t+1)-Vdx(:,t+1)').^2+(Vpy(particle,t+1)-Vdy(:,t+1)').^2+(Vpz(particle,t+1)-Vdz(:,t+1)').^2));
            probability_ij(particle,:) = (np/2+nd/2)*pi*(Dp_effective/2+Dd_effective/2)^2.*V_relative(particle,:)*10*Delta_Time_for_particles;
            minimum_distance(particle) = find(min(droplet_particle_distance(particle,:))==droplet_particle_distance(particle,:));
            dpd_min(particle) = droplet_particle_distance(particle,minimum_distance(particle));
            % A Collision happens when:
            % 1: distance dust particle to droplet < Vrelative*dtt
            % 2: probability_ij>uniform value between [0,1]
            Collision(particle,t) = dpd_min(particle)<V_relative(particle,minimum_distance(particle))...
                .*(dtt).*(probability_ij(particle,minimum_distance(particle))>rand(1,1));
        end
        mean_min_distance(t) =mean(dpd_min);    % Dummy variable to check if collisions work
        %% Collision model
        % Calculate direction vectors
        nx = -Xp(:,t)+Xd(minimum_distance,t);
        ny = -Yp(:,t)+Yd(minimum_distance,t);
        nz = -Zp(:,t)+Zd(minimum_distance,t);
        nx = nx./sqrt(nx.^2+ny.^2+nz.^2);
        ny = ny./sqrt(nx.^2+ny.^2+nz.^2);
        nz = nz./sqrt(nx.^2+ny.^2+nz.^2);
        tx = Xp(:,t)+Xd(minimum_distance,t);
        ty = Yp(:,t)+Yd(minimum_distance,t);
        tz = Zp(:,t)+Zd(minimum_distance,t); 
        tx = tx./sqrt(tx.^2+ty.^2+tz.^2);
        ty = ty./sqrt(tx.^2+ty.^2+tz.^2);
        tz = tz./sqrt(tx.^2+ty.^2+tz.^2);  
        Gbx = Vpx(:,t+1)-Vdx(minimum_distance,t+1);
        Gby = Vpy(:,t+1)-Vdy(minimum_distance,t+1);
        Gbz = Vpz(:,t+1)-Vdz(minimum_distance,t+1);
        % Re-compute the velocites after collision
        Vpx(:,t+1) = (Collision(:,t)==1).*Vpx(:,t+1)-(nx-mu_d*tx).*(nx.*Gbx+ny.*Gby+nz.*Gbz)...
            .*(1+e_r)*rho_d*Vd/(rho_p*Vp+rho_d*Vd)+(Collision(:,t)==0).*Vpx(:,t+1);
        Vpy(:,t+1) = (Collision(:,t)==1).*Vpy(:,t+1)-(ny-mu_d*ty).*(nx.*Gbx+ny.*Gby+nz.*Gbz)...
            .*(1+e_r)*rho_d*Vd/(rho_p*Vp+rho_d*Vd)+(Collision(:,t)==0).*Vpy(:,t+1);
        Vpz(:,t+1) = (Collision(:,t)==1).*Vpz(:,t+1)-(nz-mu_d*tz).*(nx.*Gbx+ny.*Gby+nz.*Gbz)...
            .*(1+e_r)*rho_d*Vd/(rho_p*Vp+rho_d*Vd)+(Collision(:,t)==0).*Vpz(:,t+1);
    else
         for particle=1:Np
            %% Landing time registration (done in this loop to reduce time
            if Zp(particle,t+1)==0 && Zp(particle,t)>0
                Landing_time(particle,1) = (t+1)*Delta_t;
                Landing_time(particle,2) = t+1;
            end
            %% Check if particle doesn't move when it touched the ground
            Xp(particle,t+1) = Xp(particle,t+1).*((Landing_time(particle,2))==1)+Xp(particle,Landing_time(particle,2)).*((Landing_time(particle,2))>1);
            Yp(particle,t+1) = Yp(particle,t+1).*((Landing_time(particle,2))==1)+Yp(particle,Landing_time(particle,2)).*((Landing_time(particle,2))>1);
            Zp(particle,t+1) = Zp(particle,t+1).*((Landing_time(particle,2))==1)+Zp(particle,Landing_time(particle,2)).*((Landing_time(particle,2))>1);
         end
    end
    % Small chance particles get stuck in rain
    Zp(:,t+1) = (Collision(:,t)==0).*Zp(:,t+1)+(Collision(:,t)==1).*(rand(Np,1)>Chance_to_get_stuck).*Zp(:,t+1);
    Time(t+1) = t*Delta_t;
    toc
end
name = ['rain_' num2str(rain_on) '_dust_particles_diameter_' num2str(Dp) '_dust_particles_number_' num2str(Np) '_rain_particles_' num2str(Nd) '_Delta_t_' num2str(Delta_t) '_Delta_t_of_particles_' num2str(Delta_Time_for_particles)];
save name
%distance = mean(mean_min_distance)
%collisions = sum(sum(Collision))