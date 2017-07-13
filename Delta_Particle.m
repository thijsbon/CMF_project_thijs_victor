%% Create variables
distance_to_sub_grid = zeros(Np,Time_steps_for_particles);
F_stokes_x_tt = zeros(Np,Time_steps_for_particles); F_stokes_x_tt(:,1) = F_stokes_x(:,t);
F_stokes_y_tt = zeros(Np,Time_steps_for_particles); F_stokes_y_tt(:,1) = F_stokes_y(:,t);
F_stokes_z_tt = zeros(Np,Time_steps_for_particles); F_stokes_z_tt(:,1) = F_stokes_z(:,t);
F_gravity_tt = zeros(Np,Time_steps_for_particles); F_gravity_tt(:,1) = F_gravity(:,t);
a_x_tt = zeros(Np,Time_steps_for_particles); a_x_tt(:,1) = a_x(:,t);
a_y_tt = zeros(Np,Time_steps_for_particles); a_y_tt(:,1) = a_y(:,t);
a_z_tt = zeros(Np,Time_steps_for_particles); a_z_tt(:,1) = a_z(:,t);
Vpx_tt = zeros(Np,Time_steps_for_particles); Vpx_tt(:,1) = Vpx(:,t);
Vpy_tt = zeros(Np,Time_steps_for_particles); Vpy_tt(:,1) = Vpy(:,t);
Vpz_tt = zeros(Np,Time_steps_for_particles); Vpz_tt(:,1) = Vpz(:,t);
Xp_tt = zeros(Np,Time_steps_for_particles); Xp_tt(:,1) = Xp(:,t);
Yp_tt = zeros(Np,Time_steps_for_particles); Yp_tt(:,1) = Yp(:,t);
Zp_tt = zeros(Np,Time_steps_for_particles); Zp_tt(:,1) = Zp(:,t);
u_prime_tt = zeros(Np,Time_steps_for_particles); u_prime_tt(:,1) = u_prime(:,t);
v_prime_tt = zeros(Np,Time_steps_for_particles); v_prime_tt(:,1) = v_prime(:,t);
w_prime_tt = zeros(Np,Time_steps_for_particles); w_prime_tt(:,1) = w_prime(:,t);
V_t = zeros(Np,Time_steps_for_particles);
eddy_life_time_tt = zeros(Np,Time_steps_for_particles); eddy_life_time_tt(:,1) = eddy_life_time(:,t);
residence_time_x = zeros(Np,Time_steps_for_particles);
residence_time_y = zeros(Np,Time_steps_for_particles);
residence_time_z = zeros(Np,Time_steps_for_particles);
a_eddy_x = zeros(Np,Time_steps_for_particles); 
a_eddy_y = zeros(Np,Time_steps_for_particles); 
a_eddy_z = zeros(Np,Time_steps_for_particles);
b_eddy = zeros(Np,Time_steps_for_particles);
k_eddy = zeros(Np,Time_steps_for_particles);
F_lift_z_tt = zeros(Np,Time_steps_for_particles);
C_L = zeros(Np,Time_steps_for_particles);
inter_particle_distance = zeros(Np,Np);
probability_ij = zeros(Np,Np);

if rain_on == 1
    ad_x_tt = zeros(Nd,Time_steps_for_particles); ad_x_tt(:,1) = ad_x(:,t);
    ad_y_tt = zeros(Nd,Time_steps_for_particles); ad_y_tt(:,1) = ad_y(:,t);
    ad_z_tt = zeros(Nd,Time_steps_for_particles); ad_z_tt(:,1) = ad_z(:,t);
    Vdx_tt = zeros(Nd,Time_steps_for_particles); Vdx_tt(:,1) = Vdx(:,t);
    Vdy_tt = zeros(Nd,Time_steps_for_particles); Vdy_tt(:,1) = Vdy(:,t);
    Vdz_tt = zeros(Nd,Time_steps_for_particles); Vdz_tt(:,1) = Vdz(:,t);
    Xd_tt = zeros(Nd,Time_steps_for_particles); Xd_tt(:,1) = Xd(:,t);
    Yd_tt = zeros(Nd,Time_steps_for_particles); Yd_tt(:,1) = Yd(:,t);
    Zd_tt = zeros(Nd,Time_steps_for_particles); Zd_tt(:,1) = Zd(:,t);
    F_stokes_xd_tt = zeros(Nd,Time_steps_for_particles); F_stokes_xd_tt(:,1) = F_stokes_xd(:,t);
    F_stokes_yd_tt = zeros(Nd,Time_steps_for_particles); F_stokes_yd_tt(:,1) = F_stokes_yd(:,t);
    F_stokes_zd_tt = zeros(Nd,Time_steps_for_particles); F_stokes_zd_tt(:,1) = F_stokes_zd(:,t);
end


%% Sub simulation
for tt = 1:Time_steps_for_particles
     %% Dust Particles
     for particle=1:Np
        %distance_to_grid(particle,t) = find(min(abs(Zp(particle,t)-zc))==abs(Zp(particle,t)-zc)).*(Zp(particle,t)>0) + (Zp(particle,t)<0);
        distance_to_sub_grid(particle,tt) = max(find(min(abs(Zp_tt(particle,tt)-zc))==abs(Zp_tt(particle,tt)-zc)));
        %distance_to_sub_grid(particle,tt) = find(min(abs(Zp_tt(particle,tt)-zc))==abs(Zp_tt(particle,tt)-zc), 1, 'first' );
        if distance_to_sub_grid(particle,tt)==Nz+2
            distance_to_sub_grid(particle,tt)=Nz+1;
        end
        if distance_to_sub_grid(particle,tt)==1
            distance_to_sub_grid(particle,tt)=2;
        end
     end
     
     %% Eddy simulation
     if Euler_Lagrangian_Eddy == 1
        %% Langevin model
        V_t(:,tt) = (Zp_tt(:,tt)<500).*l_effective(distance_to_sub_grid(:,tt),t).*(u(distance_to_sub_grid(:,tt)+1,t)-u(distance_to_sub_grid(:,tt)-1,t))./zc(distance_to_sub_grid(:,tt))' +l_effective(end,t).*(u(end-1,t)-u(end-2,t))./zc(end)'.*(Zp_tt(:,tt)>500);
        eddy_life_time_tt(:,tt)  = c_t * l_effective(distance_to_sub_grid(:,tt))./V_t(:,tt);
        residence_time_x(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpx_tt(:,tt)-u(distance_to_sub_grid(:,tt),t)));
        residence_time_y(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpy_tt(:,tt)-u(distance_to_sub_grid(:,tt),t)));
        residence_time_z(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpz_tt(:,tt)-u(distance_to_sub_grid(:,tt),t)));
        Ti_star(:,tt) = l_effective(distance_to_sub_grid(:,tt),t)./V_t(:,tt);
        a_eddy_x(:,tt) = exp(-Delta_Time_for_particles./Ti_star(:,tt));
        a_eddy_y(:,tt) = exp(-Delta_Time_for_particles./Ti_star(:,tt));
        a_eddy_z(:,tt) = exp(-Delta_Time_for_particles./Ti_star(:,tt));
        %k_eddy(:,tt) = sqrt(u_prime_tt(:,tt).^2+v_prime_tt(:,tt).^2+w_prime_tt(:,tt).^2);
        k_eddy(:,tt) = V_t(:,tt)/c_t;
        b_eddy(:,tt) = k_eddy(:,tt).*sqrt(1-mean([a_eddy_x(:,tt);a_eddy_y(:,tt);a_eddy_z(:,tt)]).^2);
        u_prime_tt(:,tt+1) = a_eddy_x(:,tt).*u_prime_tt(:,tt)+b_eddy(:,tt).*randn(Np,1);
        v_prime_tt(:,tt+1) = a_eddy_y(:,tt).*v_prime_tt(:,tt)+b_eddy(:,tt).*randn(Np,1);
        w_prime_tt(:,tt+1) = a_eddy_z(:,tt).*w_prime_tt(:,tt)+b_eddy(:,tt).*randn(Np,1);
     elseif Euler_Lagrangian_Eddy == 2
        %% Descrete eddy model
        V_t(:,tt) = l_effective(distance_to_sub_grid(:,tt),t).*(u(distance_to_sub_grid(:,tt)+1,t)-u(distance_to_sub_grid(:,tt)-1,t))./zc(distance_to_sub_grid(:,tt))';
        eddy_life_time_tt(:,tt)  = c_t * l_effective(distance_to_sub_grid(:,tt),t)./V_t(:,tt);
        residence_time_x(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpx_tt(:,tt)-u_prime_tt(:,tt)-u(distance_to_sub_grid(:,tt),t)));
        residence_time_y(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpy_tt(:,tt)-v_prime_tt(:,tt)));
        residence_time_z(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpz_tt(:,tt)-w_prime_tt(:,tt)));
        %residence_time_y(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpy_tt(:,t)-u(distance_to_sub_grid(:,tt),t)));
        %residence_time_z(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpz_tt(:,t)-u(distance_to_sub_grid(:,tt),t)));
        T_I_x(:,tt) = min(eddy_life_time_tt(:,tt),residence_time_x(:,tt));
        T_I_y(:,tt) = min(eddy_life_time_tt(:,tt),residence_time_y(:,tt));
        T_I_z(:,tt) = min(0.1,min(eddy_life_time_tt(:,tt),residence_time_z(:,tt)));
        %k_eddy(:,tt) = sqrt(u_prime_tt(:,tt).^2+v_prime_tt(:,tt).^2+w_prime_tt(:,tt).^2);
        k_eddy(:,tt) = 10*V_t(:,tt)/c_t;
        u_prime_tt(:,tt+1) = u_prime_tt(:,tt).*(t*Delta_t<Eddy_start_x+T_I_x(:,tt)) + k_eddy(:,tt).*randn(Np,1).*(t*Delta_t>Eddy_start_x+T_I_x(:,tt));
        v_prime_tt(:,tt+1) = v_prime_tt(:,tt).*(t*Delta_t<Eddy_start_y+T_I_y(:,tt)) + k_eddy(:,tt).*randn(Np,1).*(t*Delta_t>Eddy_start_y+T_I_y(:,tt));
        w_prime_tt(:,tt+1) = w_prime_tt(:,tt).*(t*Delta_t<Eddy_start_z+T_I_z(:,tt)) + k_eddy(:,tt).*randn(Np,1).*(t*Delta_t>Eddy_start_z+T_I_z(:,tt));
        % New Eddy start times
        Eddy_start_x = (tt*Delta_Time_for_particles>Eddy_start_x+T_I_x(:,tt)).*tt*Delta_Time_for_particles;
        Eddy_start_y = (tt*Delta_Time_for_particles>Eddy_start_y+T_I_y(:,tt)).*tt*Delta_Time_for_particles;
        Eddy_start_z = (tt*Delta_Time_for_particles>Eddy_start_z+T_I_z(:,tt)).*tt*Delta_Time_for_particles;
     end
     % Check if u,v,w prime stay reasonable 
     u_prime_tt(:,tt+1) = (abs(u_prime_tt(:,tt+1))<6.5).*u_prime_tt(:,tt+1)+(u_prime_tt(:,tt+1)>6.5)*6.5+(u_prime_tt(:,tt+1)<-15)*-6.5;
     v_prime_tt(:,tt+1) = (abs(v_prime_tt(:,tt+1))<6.5).*v_prime_tt(:,tt+1)+(v_prime_tt(:,tt+1)>6.5)*6.5+(v_prime_tt(:,tt+1)<-15)*-6.5;
     w_prime_tt(:,tt+1) = (abs(w_prime_tt(:,tt+1))<6.5).*w_prime_tt(:,tt+1)+(w_prime_tt(:,tt+1)>6.5)*6.5+(w_prime_tt(:,tt+1)<-15)*-6.5;     
     
     %% Classical Mechanics
     % Calculate forces
     F_stokes_x_tt(:,tt) = C_stokes * (u(distance_to_sub_grid(:,tt),t)+u_prime_tt(:,tt+1)-Vpx_tt(:,tt));
     F_stokes_y_tt(:,tt) = C_stokes * (v_prime_tt(:,tt+1)-Vpy_tt(:,tt));
     F_stokes_z_tt(:,tt) = C_stokes * (w_prime_tt(:,tt+1)-Vpz_tt(:,tt));
     %F_gravity_tt(:,tt) = -g*rho_p*Vp;
     F_stokes_x_tt(:,tt) = (Zp_tt(:,tt)>0).*F_stokes_x_tt(:,tt);
     F_stokes_y_tt(:,tt) = (Zp_tt(:,tt)>0).*F_stokes_y_tt(:,tt);
     F_stokes_z_tt(:,tt) = (Zp_tt(:,tt)>0).*F_stokes_z_tt(:,tt);
     %F_gravity_tt(:,tt) = (Zp_tt(:,tt)>0).*F_gravity_tt(:,t);
     % Calculate accelerations 
     a_x_tt(:,tt) = F_stokes_x_tt(:,tt)/(Vp*rho_p);
     a_y_tt(:,tt) = F_stokes_y_tt(:,tt)/(Vp*rho_p);
     %a_z_tt(:,tt) = (F_stokes_z_tt(:,tt))/(Vp*rho_p)-g*Dp/2*rho_p*0.1;
     a_z_tt(:,tt) = (F_stokes_z_tt(:,tt))/(Vp*rho_p)-g;
     % Calculate velocities OLD
     Vpx_tt(:,tt+1) = Vpx_tt(:,tt) + a_x_tt(:,tt)*Delta_Time_for_particles;
     Vpy_tt(:,tt+1) = Vpy_tt(:,tt) + a_y_tt(:,tt)*Delta_Time_for_particles;
     Vpz_tt(:,tt+1) = Vpz_tt(:,tt) + a_z_tt(:,tt)*Delta_Time_for_particles;
     
     % Calculate new positions
     Xp_tt(:,tt+1) = Xp_tt(:,tt) + Vpx_tt(:,tt)*Delta_Time_for_particles + 0.5*a_x_tt(:,tt)*Delta_Time_for_particles^2;
     Yp_tt(:,tt+1) = Yp_tt(:,tt) + Vpy_tt(:,tt)*Delta_Time_for_particles + 0.5*a_y_tt(:,tt)*Delta_Time_for_particles^2;
     Zp_tt(:,tt+1) = Zp_tt(:,tt) + Vpz_tt(:,tt)*Delta_Time_for_particles + 0.5*a_z_tt(:,tt)*Delta_Time_for_particles^2;
     % Check if particle does not go below ground
     Zp_tt(:,tt+1) = (~isnan(Zp_tt(:,tt+1))).*(Zp_tt(:,tt+1)>0).*Zp_tt(:,tt+1);
     Vpx_tt(:,tt+1) = (Zp_tt(:,tt+1)>0).*Vpx_tt(:,tt+1);
     Vpy_tt(:,tt+1) = (Zp_tt(:,tt+1)>0).*Vpy_tt(:,tt+1);
     Vpz_tt(:,tt+1) = (Zp_tt(:,tt+1)>0).*Vpz_tt(:,tt+1);
     
     % Check if velocities stay realisitc
     %Vpx_tt(:,tt+1) = (Vpx_tt(:,tt+1)<20).*Vpx_tt(:,tt+1)+(Vpx_tt(:,tt+1)>20)*20;
     %Vpy_tt(:,tt+1) = (Vpy_tt(:,tt+1)<15).*Vpy_tt(:,tt+1)+(Vpy_tt(:,tt+1)>15)*15;
     %Vpz_tt(:,tt+1) = (Vpz_tt(:,tt+1)<20).*Vpz_tt(:,tt+1)+(Vpz_tt(:,tt+1)>20)*20;
     
     %% Inter particle collisions NIET AF!
%      for particle=1:Np
%          inter_particle_distance(particle,:) = sqrt((Xp_tt(:,tt)-Xp_tt(particle,tt)).^2+(Yp_tt(:,tt)-Yp_tt(particle,tt)).^2+(Zp_tt(:,tt)-Zp_tt(particle,tt)).^2)';
%          if inter_particle_distance(particle,:)'<sqrt((Vpx_tt(:,tt)-Vpx_tt(particle,tt)).^2+(Vpy_tt(:,tt)-Vpy_tt(particle,tt)).^2+(Vpz_tt(:,tt)-Vpz_tt(particle,tt)).^2)*Delta_Time_for_particles
%             probability_ij(particle,:) = np*pi*Dp_effective^2.*sqrt((Vpx_tt(:,tt)-Vpx_tt(particle,tt)).^2+(Vpy_tt(:,tt)-Vpy_tt(particle,tt)).^2+(Vpz_tt(:,tt)-Vpz_tt(particle,tt)).^2)*Delta_Time_for_particles;
%          end
%      end

     %% Droplets
     if rain_on == 1
         F_stokes_xd_tt(:,tt) = C_stokes_rain * (max(u(:,t))-Vdx_tt(:,tt));
         F_stokes_yd_tt(:,tt) = C_stokes_rain * (-Vdy_tt(:,tt));
         F_stokes_zd_tt(:,tt) = C_stokes_rain * (-Vdz_tt(:,tt));
         ad_x_tt(:,tt) = F_stokes_xd_tt(:,tt)/(rho_d*Vd);
         ad_y_tt(:,tt) = F_stokes_yd_tt(:,tt)/(rho_d*Vd);
         ad_z_tt(:,tt) = F_stokes_zd_tt(:,tt)/(rho_d*Vd)-g;
         Vdx_tt(:,tt+1) = Vdx_tt(:,tt)+ad_x_tt(:,tt)*Delta_Time_for_particles;
         Vdy_tt(:,tt+1) = Vdy_tt(:,tt)+ad_y_tt(:,tt)*Delta_Time_for_particles;
         Vdz_tt(:,tt+1) = Vdz_tt(:,tt)+ad_z_tt(:,tt)*Delta_Time_for_particles;
    
         Xd_tt(:,tt+1) = Xd_tt(:,tt)+Vdx_tt(:,tt)*Delta_Time_for_particles+0.5*ad_x_tt(:,tt)*Delta_Time_for_particles^2;
         Yd_tt(:,tt+1) = Yd_tt(:,tt)+Vdx_tt(:,tt)*Delta_Time_for_particles+0.5*ad_y_tt(:,tt)*Delta_Time_for_particles^2;
         Zd_tt(:,tt+1) = Zd_tt(:,tt)+Vdx_tt(:,tt)*Delta_Time_for_particles+0.5*ad_z_tt(:,tt)*Delta_Time_for_particles^2;
         
         % Go back to 1 m^3
         Xd_tt(:,tt+1) = Xd_tt(:,tt+1)-floor(Xd_tt(:,tt+1));
         Yd_tt(:,tt+1) = Yd_tt(:,tt+1)-floor(Yd_tt(:,tt+1));
         Zd_tt(:,tt+1) = Zd_tt(:,tt+1)-floor(Zd_tt(:,tt+1));
     end
end

if speed_up == 1
    %% Statistics
    % Dust Particles
    F_stokes_x_statistics = [mean(F_stokes_x_tt')' std(F_stokes_x_tt')'];
    F_stokes_y_statistics = [mean(F_stokes_y_tt')' std(F_stokes_y_tt')'];
    F_stokes_z_statistics = [mean(F_stokes_z_tt')' std(F_stokes_z_tt')'];
    a_x_statistics = [mean(a_x_tt')' std(a_x_tt')'];
    a_y_statistics = [mean(a_y_tt')' std(a_y_tt')'];
    a_z_statistics = [mean(a_z_tt')' std(a_z_tt')'];
    Vpx_statistics = [mean(Vpx_tt')' std(Vpx_tt')'];
    Vpy_statistics = [mean(Vpy_tt')' std(Vpy_tt')'];
    Vpz_statistics = [mean(Vpz_tt')' std(Vpz_tt')'];
    Xp_statistics = [mean(Xp_tt')' std(Xp_tt')'];
    Yp_statistics = [mean(Yp_tt')' std(Yp_tt')'];
    Zp_statistics = [mean(Zp_tt')' std(Zp_tt')'];
    u_prime_statistics = [mean(u_prime_tt')' std(u_prime_tt')'];
    v_prime_statistics = [mean(v_prime_tt')' std(v_prime_tt')'];
    w_prime_statistics = [mean(w_prime_tt')' std(w_prime_tt')'];
    eddy_life_time_statistics = [mean(eddy_life_time_tt')' std(eddy_life_time_tt')'];
    if rain_on ==1
    % Droplet
    ad_x_statistics = [mean(ad_x_tt')' std(ad_x_tt')'];
    ad_y_statistics = [mean(ad_y_tt')' std(ad_y_tt')'];
    ad_z_statistics = [mean(ad_z_tt')' std(ad_z_tt')'];
    Vdx_statistics = [mean(Vdx_tt')' std(Vdx_tt')'];
    Vdy_statistics = [mean(Vdy_tt')' std(Vdy_tt')'];
    Vdz_statistics = [mean(Vdz_tt')' std(Vdz_tt')'];
    Xd_statistics = [mean(Xd_tt')' std(Xd_tt')'];
    Yd_statistics = [mean(Yd_tt')' std(Yd_tt')'];
    Zd_statistics = [mean(Zd_tt')' std(Zd_tt')'];
    F_stokes_xd_statistics = [mean(F_stokes_xd_tt')' std(F_stokes_xd_tt')'];
    F_stokes_yd_statistics = [mean(F_stokes_yd_tt')' std(F_stokes_yd_tt')'];
    F_stokes_zd_statistics = [mean(F_stokes_zd_tt')' std(F_stokes_zd_tt')'];
    end
    %% New end values
    % Dust particles
    F_stokes_x_tt(:,end) = F_stokes_x_statistics(:,1)+F_stokes_x_statistics(:,2).*randn(Np,1);
    F_stokes_y_tt(:,end) = F_stokes_y_statistics(:,1)+F_stokes_y_statistics(:,2).*randn(Np,1);
    F_stokes_z_tt(:,end) = F_stokes_z_statistics(:,1)+F_stokes_z_statistics(:,2).*randn(Np,1);
    a_x_tt(:,end) = a_x_statistics(:,1)+a_x_statistics(:,2).*randn(Np,1);
    a_y_tt(:,end) = a_y_statistics(:,1)+a_y_statistics(:,2).*randn(Np,1);
    a_z_tt(:,end) = a_z_statistics(:,1)+a_z_statistics(:,2).*randn(Np,1);
    u_prime_tt(:,end) = u_prime_statistics(:,1)+u_prime_statistics(:,2).*randn(Np,1);
    v_prime_tt(:,end) = v_prime_statistics(:,1)+v_prime_statistics(:,2).*randn(Np,1);
    w_prime_tt(:,end) = w_prime_statistics(:,1)+w_prime_statistics(:,2).*randn(Np,1);
    if speed_up_method == 1
        Vpx_tt(:,end) = Vpx_statistics(:,1)+Vpx_statistics(:,2).*randn(Np,1);%+a_x_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles;
        Vpy_tt(:,end) = Vpy_statistics(:,1)+Vpy_statistics(:,2).*randn(Np,1);%+a_y_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles;
        Vpz_tt(:,end) = Vpz_statistics(:,1)+Vpz_statistics(:,2).*randn(Np,1);%+a_z_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles;
        Xp_tt(:,end) = Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1);
        Yp_tt(:,end) = Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1);
        Zp_tt(:,end) = Zp_statistics(:,1)+Zp_statistics(:,2).*randn(Np,1);
    elseif speed_up_method == 2
        Vpx_tt(:,end) = Vpx(:,end-1)+a_x_tt(:,end).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles;
        Vpy_tt(:,end) = Vpy(:,end-1)+a_y_tt(:,end).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles;
        Vpz_tt(:,end) = Vpz(:,end-1)+a_z_tt(:,end).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles;
        if sub_method_2 == 1
            Xp_tt(:,end) = Xp_tt(:,end-1)+Vpx(:,end-1).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles+0.5*a_x_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles^2;
            Yp_tt(:,end) = Yp_tt(:,end-1)+Vpy(:,end-1).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles+0.5*a_y_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles^2;
            Zp_tt(:,end) = Zp_tt(:,end-1)+Vpz(:,end-1).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles+0.5*a_z_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles^2;
        elseif sub_method_2 == 2
            Xp_tt(:,end) = Xp_tt(:,end-1)+Vpx(:,end).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles+0.5*a_x_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles^2;
            Yp_tt(:,end) = Yp_tt(:,end-1)+Vpy(:,end).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles+0.5*a_y_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles^2;
            Zp_tt(:,end) = Zp_tt(:,end-1)+Vpz(:,end).*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles+0.5*a_z_tt(:,end)*(Original_Time_steps_for_particles-Time_steps_for_particles)*Delta_Time_for_particles^2;
        elseif sub_method_2 == 3
            Xp_tt(:,end) = Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1)+Vpx_tt(:,end-1)*Delta_t;
            Yp_tt(:,end) = Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1)+Vpy_tt(:,end-1)*Delta_t;
            Zp_tt(:,end) = Zp_statistics(:,1)+Zp_statistics(:,2).*randn(Np,1)+Vpz_tt(:,end-1)*Delta_t;
        end
    elseif speed_up_method == 3
        Vpx_tt(:,end) = Vpx_statistics(:,1)+Vpx_statistics(:,2).*randn(Np,1);
        Vpy_tt(:,end) = Vpy_statistics(:,1)+Vpy_statistics(:,2).*randn(Np,1);
        Vpz_tt(:,end) = Vpz_statistics(:,1)+Vpz_statistics(:,2).*randn(Np,1);
        Xp_tt(:,end) = ((Xp_tt(:,end)-Xp_tt(:,1))<1*Xp_statistics(:,2))...
            .*((Xp_tt(:,end)-Xp_tt(:,1))/(Time_steps_for_particles*Delta_Time_for_particles-Delta_Time_for_particles)...
            .*(Delta_t-Delta_Time_for_particles)+Xp_tt(:,1))...
            +((Xp_tt(:,end)-Xp_tt(:,1))>2*Xp_statistics(:,1)).*(Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1));
        Yp_tt(:,end) = ((Yp_tt(:,end)-Yp_tt(:,1))<1*Yp_statistics(:,2))...
            .*((Yp_tt(:,end)-Yp_tt(:,1))/(Time_steps_for_particles*Delta_Time_for_particles-Delta_Time_for_particles)...
            .*(Delta_t-Delta_Time_for_particles)+Yp_tt(:,1))...
            +((Yp_tt(:,end)-Yp_tt(:,1))>2*Yp_statistics(:,1)).*(Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1));
        Zp_tt(:,end) = ((Zp_tt(:,end)-Zp_tt(:,1))<1*Zp_statistics(:,2))...
            .*((Zp_tt(:,end)-Zp_tt(:,1))/(Time_steps_for_particles*Delta_Time_for_particles-Delta_Time_for_particles)...
            .*(Delta_t-Delta_Time_for_particles)+Zp_tt(:,1))...
            +((Zp_tt(:,end)-Zp_tt(:,1))>2*Zp_statistics(:,1)).*(Zp_statistics(:,1)+Zp_statistics(:,2).*randn(Np,1));
    elseif speed_up_method == 4
        Vpx_tt(:,end) = Vpx_statistics(:,1)+Vpx_statistics(:,2).*randn(Np,1);
        Vpy_tt(:,end) = Vpy_statistics(:,1)+Vpy_statistics(:,2).*randn(Np,1);
        Vpz_tt(:,end) = Vpz_statistics(:,1)+Vpz_statistics(:,2).*randn(Np,1);
        X_method_4_1 = ((Xp_tt(:,end)-Xp_tt(:,1))<1*Xp_statistics(:,2))...
            .*((Xp_tt(:,end)-Xp_tt(:,1))/(Time_steps_for_particles*Delta_Time_for_particles-Delta_Time_for_particles)...
            .*(Delta_t-Delta_Time_for_particles)+Xp_tt(:,1))...
            +((Xp_tt(:,end)-Xp_tt(:,1))>2*Xp_statistics(:,1)).*(Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1));
        X_method_4_2 = Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1);
        Y_method_4_1 = ((Yp_tt(:,end)-Yp_tt(:,1))<1*Yp_statistics(:,2))...
            .*((Yp_tt(:,end)-Yp_tt(:,1))/(Time_steps_for_particles*Delta_Time_for_particles-Delta_Time_for_particles)...
            .*(Delta_t-Delta_Time_for_particles)+Yp_tt(:,1))...
            +((Yp_tt(:,end)-Yp_tt(:,1))>2*Yp_statistics(:,1)).*(Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1));
        Y_method_4_2 = Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1);
        Z_method_4_1 = ((Zp_tt(:,end)-Zp_tt(:,1))<1*Zp_statistics(:,2))...
            .*((Zp_tt(:,end)-Zp_tt(:,1))/(Time_steps_for_particles*Delta_Time_for_particles-Delta_Time_for_particles)...
            .*(Delta_t-Delta_Time_for_particles)+Zp_tt(:,1))...
            +((Zp_tt(:,end)-Zp_tt(:,1))>2*Zp_statistics(:,1)).*(Zp_statistics(:,1)+Zp_statistics(:,2).*randn(Np,1));
        Z_method_4_2 = Zp_statistics(:,1)+Zp_statistics(:,2).*randn(Np,1);
        Xp_tt(:,end) = (X_method_4_1>0).*(X_method_4_2>0).*mean([X_method_4_1'; X_method_4_2'])';
        Yp_tt(:,end) = (Y_method_4_1>0).*(Y_method_4_2>0).*mean([Y_method_4_1'; Y_method_4_2'])';
        Zp_tt(:,end) = (Z_method_4_1>0).*(Z_method_4_2>0).*mean([Z_method_4_1'; Z_method_4_2'])';
    elseif speed_up_method == 5
        Vpx_tt(:,end) = Vpx_statistics(:,1)+Vpx_statistics(:,2).*randn(Np,1);
        Vpy_tt(:,end) = Vpy_statistics(:,1)+Vpy_statistics(:,2).*randn(Np,1);
        Vpz_tt(:,end) = Vpz_statistics(:,1)+Vpz_statistics(:,2).*randn(Np,1);
        Xp_tt(:,end) = Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1)+Vpx_statistics(:,1)*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
        Yp_tt(:,end) = Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1)+Vpy_statistics(:,1)*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
        Zp_tt(:,end) = Zp_tt(:,end)-2/9*(rho_p-rho)/mu*g*(Dp/2)^2*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
    elseif speed_up_method == 6
        Vpx_tt(:,end) = Vpx_statistics(:,1)+Vpx_statistics(:,2).*randn(Np,1);
        Vpy_tt(:,end) = Vpy_statistics(:,1)+Vpy_statistics(:,2).*randn(Np,1);
        Vpz_tt(:,end) = Vpz_statistics(:,1)+Vpz_statistics(:,2).*randn(Np,1);
        Xp_tt(:,end) = Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1)+Vpx_statistics(:,1)*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
        Yp_tt(:,end) = Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1)+Vpy_statistics(:,1)*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
        Zp_tt(:,end) = Zp_tt(:,end)+mod_6*Zp_statistics(:,2).*randn(Np,1)-2/9*(rho_p-rho)/mu*g*(Dp/2)^2*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
    elseif speed_up_method == 7
        Vpx_tt(:,end) = Vpx_statistics(:,1)+Vpx_statistics(:,2).*randn(Np,1);
        Vpy_tt(:,end) = Vpy_statistics(:,1)+Vpy_statistics(:,2).*randn(Np,1);
        Vpz_tt(:,end) = Vpz_statistics(:,1)+Vpz_statistics(:,2).*randn(Np,1);
        Xp_tt(:,end) = Xp_statistics(:,1)+Xp_statistics(:,2).*randn(Np,1)+add_percentage*Vpx_statistics(:,1)*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
        Yp_tt(:,end) = Yp_statistics(:,1)+Yp_statistics(:,2).*randn(Np,1)+Vpy_statistics(:,1)*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
        Zp_tt(:,end) = Zp_tt(:,end)-2/9*(rho_p-rho)/mu*g*(Dp/2)^2*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);

    end
    % Check if particle does not go below ground
    Zp_tt(:,end) = (~isnan(Zp_tt(:,end))).*(Zp_tt(:,end)>0).*Zp_tt(:,end);
    Xp_tt(:,end) = (~isnan(Zp_tt(:,end))).*(Zp_tt(:,end)>0).*Xp_tt(:,end)+(~isnan(Zp_tt(:,end))).*(Zp_tt(:,end)==0).*Xp_tt(:,tt-1);
    Yp_tt(:,end) = (~isnan(Zp_tt(:,end))).*(Zp_tt(:,end)>0).*Yp_tt(:,end)+(~isnan(Zp_tt(:,end))).*(Zp_tt(:,end)==0).*Yp_tt(:,tt-1);
    Vpx_tt(:,end) = (Zp_tt(:,end)>0).*Vpx_tt(:,end);
    Vpy_tt(:,end) = (Zp_tt(:,end)>0).*Vpy_tt(:,end);
    Vpz_tt(:,end) = (Zp_tt(:,end)>0).*Vpz_tt(:,end);
    
    if rain_on==1
    % Droplet(Not really necessary)
    Droplet_necessary =1;
    if Droplet_necessary == 1
        ad_x_tt(:,end) = ad_x_statistics(:,1)+ad_x_statistics(:,2).*randn(Nd,1);
        ad_y_tt(:,end) = ad_y_statistics(:,1)+ad_y_statistics(:,2).*randn(Nd,1);
        ad_z_tt(:,end) = ad_z_statistics(:,1)+ad_z_statistics(:,2).*randn(Nd,1);
        Vdx_tt(:,end) = Vdx_statistics(:,1)+Vdx_statistics(:,2).*randn(Nd,1);
        Vdy_tt(:,end) = Vdy_statistics(:,1)+Vdy_statistics(:,2).*randn(Nd,1)+0.01*randn(Nd,1);
        Vdz_tt(:,end) = Vdz_statistics(:,1)+Vdz_statistics(:,2).*randn(Nd,1);
        %Vdz_tt(:,end) = -2/9*(rho_d-rho)/mu*g*(Dd/2)^2+Vdz_statistics(:,2).*randn(Nd,1);
        Xd_tt(:,end) = Xd_statistics(:,1)+Xd_statistics(:,2).*randn(Nd,1);
        Yd_tt(:,end) = Yd_statistics(:,1)+Yd_statistics(:,2).*randn(Nd,1);
        Zd_tt(:,end) = Zd_tt(:,end)-2/9*(rho_d-rho)/mu*g*(Dd/2)^2*(Delta_t-Time_steps_for_particles*Delta_Time_for_particles);
        F_stokes_xd_tt(:,end) = F_stokes_xd_statistics(:,1)+F_stokes_xd_statistics(:,2).*randn(Nd,1);
        F_stokes_yd_tt(:,end) = F_stokes_yd_statistics(:,1)+F_stokes_yd_statistics(:,2).*randn(Nd,1);
        F_stokes_zd_tt(:,end) = F_stokes_zd_statistics(:,1)+F_stokes_zd_statistics(:,2).*randn(Nd,1);
        % Go back to 1 m^3
        Xd_tt(:,end) = Xd_tt(:,end)-floor(Xd_tt(:,end));
        Yd_tt(:,end) = Yd_tt(:,end)-floor(Yd_tt(:,end));
        Zd_tt(:,end) = Zd_tt(:,end)-floor(Zd_tt(:,end));
    end
    end
end