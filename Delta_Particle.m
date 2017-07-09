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

%% Sub simulation
for tt = 1:Time_steps_for_particles
     for particle=1:Np
        %distance_to_grid(particle,t) = find(min(abs(Zp(particle,t)-zc))==abs(Zp(particle,t)-zc)).*(Zp(particle,t)>0) + (Zp(particle,t)<0);
        distance_to_sub_grid(particle,tt) = max(find(min(abs(Zp_tt(particle,tt)-zc))==abs(Zp_tt(particle,tt)-zc)));
        %distance_to_sub_grid(particle,tt) = find(min(abs(Zp_tt(particle,tt)-zc))==abs(Zp_tt(particle,tt)-zc), 1, 'first' );
        if distance_to_sub_grid(particle,tt)==Nz+2
            distance_to_sub_grid(particle,tt)=Nz+1;
        end
     end
     
     %% Eddy simulation
     if Euler_Lagrangian_Eddy == 1
        %% Langevin model
        V_t(:,tt) = (Zp_tt(:,tt)<1000).*l_effective(distance_to_sub_grid(:,tt),t).*(u(distance_to_sub_grid(:,tt)+1,t)-u(distance_to_sub_grid(:,tt)-1,t))./zc(distance_to_sub_grid(:,tt))' + 1e5*(Zp_tt(:,tt)>1000);
        eddy_life_time_tt(:,tt)  = c_t * l_effective(distance_to_sub_grid(:,tt))./V_t(:,tt);
        residence_time_x(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpx_tt(:,tt)-u(distance_to_sub_grid(:,tt),t)));
        residence_time_y(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpy_tt(:,tt)-u(distance_to_sub_grid(:,tt),t)));
        residence_time_z(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpz_tt(:,tt)-u(distance_to_sub_grid(:,tt),t)));
        Ti_star(:,tt) = l_effective(distance_to_sub_grid(:,tt),t)./V_t(:,tt);
        a_eddy_x(:,tt) = exp(-Delta_Time_for_particles./Ti_star(:,tt));
        a_eddy_y(:,tt) = exp(-Delta_Time_for_particles./Ti_star(:,tt));
        a_eddy_z(:,tt) = exp(-Delta_Time_for_particles./Ti_star(:,tt));
        k_eddy(:,tt) = sqrt(u_prime_tt(:,tt).^2+v_prime_tt(:,tt).^2+w_prime_tt(:,tt).^2);
        b_eddy(:,tt) = k_eddy(:,tt).*sqrt(1-mean([a_eddy_x(:,tt);a_eddy_y(:,tt);a_eddy_z(:,tt)]).^2);
        u_prime_tt(:,tt+1) = a_eddy_x(:,tt).*u_prime_tt(:,tt)+b_eddy(:,tt).*randn(Np,1);
        v_prime_tt(:,tt+1) = a_eddy_y(:,tt).*v_prime_tt(:,tt)+b_eddy(:,tt).*randn(Np,1);
        w_prime_tt(:,tt+1) = a_eddy_z(:,tt).*w_prime_tt(:,tt)+b_eddy(:,tt).*randn(Np,1);
     elseif Euler_Lagrangian_Eddy == 2
        %% Descrete eddy model
        V_t(:,tt) = l_effective(distance_to_sub_grid(:,tt),t).*(u(distance_to_sub_grid(:,tt)+1,t)-u(distance_to_sub_grid(:,tt)-1,t))./zc(distance_to_sub_grid(:,tt))';
        eddy_life_time_tt(:,tt)  = c_t * l_effective(distance_to_sub_grid(:,tt),t)./V_t(:,tt);
        residence_time_x(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpx_tt(:,t)-u(distance_to_sub_grid(:,tt),t)));
        residence_time_y(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpy_tt(:,t)-u(distance_to_sub_grid(:,tt),t)));
        residence_time_z(:,tt)  = l_effective(distance_to_sub_grid(:,tt),t)./(abs(Vpz_tt(:,t)-u(distance_to_sub_grid(:,tt),t)));
        T_I_x(:,tt) = min(eddy_life_time_tt(:,tt),residence_time_x(:,tt));
        T_I_y(:,tt) = min(eddy_life_time_tt(:,tt),residence_time_y(:,tt));
        T_I_z(:,tt) = min(0.1,min(eddy_life_time_tt(:,tt),residence_time_z(:,tt)));
        u_prime_tt(:,tt+1) = u_prime_tt(:,tt).*(t*Delta_t<Eddy_start_x+T_I_x(:,tt)) + randn(Np,1).*(t*Delta_t>Eddy_start_x+T_I_x(:,tt));
        v_prime_tt(:,tt+1) = v_prime_tt(:,tt).*(t*Delta_t<Eddy_start_y+T_I_y(:,tt)) + randn(Np,1).*(t*Delta_t>Eddy_start_y+T_I_y(:,tt));
        w_prime_tt(:,tt+1) = w_prime_tt(:,tt).*(t*Delta_t<Eddy_start_z+T_I_z(:,tt)) + (1.5*randn(Np,1)+upper).*(t*Delta_t>Eddy_start_z+T_I_z(:,tt));
        % New Eddy start times
        Eddy_start_x = (tt*Delta_Time_for_particles>Eddy_start_x+T_I_x(:,tt)).*tt*Delta_Time_for_particles;
        Eddy_start_y = (tt*Delta_Time_for_particles>Eddy_start_y+T_I_y(:,tt)).*tt*Delta_Time_for_particles;
        Eddy_start_z = (tt*Delta_Time_for_particles>Eddy_start_z+T_I_z(:,tt)).*tt*Delta_Time_for_particles;
     end
     % Check if u,v,w prime stay reasonable 
     u_prime_tt(:,tt+1) = (u_prime_tt(:,tt+1)<15).*u_prime_tt(:,tt+1)+(u_prime_tt(:,tt+1)>15)*15;
     v_prime_tt(:,tt+1) = (v_prime_tt(:,tt+1)<15).*v_prime_tt(:,tt+1)+(v_prime_tt(:,tt+1)>15)*15;
     w_prime_tt(:,tt+1) = (w_prime_tt(:,tt+1)<15).*w_prime_tt(:,tt+1)+(w_prime_tt(:,tt+1)>15)*15;     
     
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
     a_z_tt(:,tt) = (F_stokes_z_tt(:,tt))/(Vp*rho_p)-g*Dp/2*rho_p*0.1;
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
     Vpz_tt(:,tt+1) = (Zp_tt(:,tt+1)>0).*Vpz_tt(:,tt+1);
     
     % Check if velocities stay realisitc
     Vpx_tt(:,tt+1) = (Vpx_tt(:,tt+1)<20).*Vpx_tt(:,tt+1)+(Vpx_tt(:,tt+1)>20)*20;
     Vpy_tt(:,tt+1) = (Vpy_tt(:,tt+1)<15).*Vpy_tt(:,tt+1)+(Vpy_tt(:,tt+1)>15)*15;
     Vpz_tt(:,tt+1) = (Vpz_tt(:,tt+1)<20).*Vpz_tt(:,tt+1)+(Vpz_tt(:,tt+1)>20)*20;
     
     %% Inter particle collisions NIET AF!
     for particle=1:Np
        probability_ij_x(particle,:) = np*pi*Dp_effective^2*(Vpx_tt(:,tt)-Vpx(particle,tt))*Delta_Time_for_particles;
        probability_ij_y(particle,:) = np*pi*Dp_effective^2*(Vpy_tt(:,tt)-Vpy(particle,tt))*Delta_Time_for_particles;
        probability_ij_z(particle,:) = np*pi*Dp_effective^2*(Vpz_tt(:,tt)-Vpz(particle,tt))*Delta_Time_for_particles;
        probability_ij(particle,:) = probability_ij_x.*probability_ij_y.*probability_ij_z;
     end
     
end