iter = 0;
residue = 10*min_residue;
Qnew = 0;
u_change =1;
while iter<max_iter && residue>min_residue
    %% velocity change (for residue determination)
    
    
    %% dpdx change (in case of the prescribed flow rate)
    %dpdx = dpdx+(max_iter-iter)*mean(nu)*(Qnew-Q)/H^3*(prescribeswitch == 1);
    dpdx = dpdx + (Qnew-Q)/Q*0.01*(prescribeswitch == 1);
    
    %%    
    %tau_wall_D = abs(u(2)-u(1))/dzc(1); %Tau wall down
    %tau_wall_U = abs(u(end)-u(end-1))/dzc(end); %Tau wall up
    %u_star_U = sqrt(tau_wall_U/rho);
    %u_star_D = sqrt(tau_wall_D/rho);
    %y_star_D = zc(2:end-1).*u_star_D./nu_c;
    %y_star_U = (H-zc(2:end-1)* (bcswitch==0) ).*u_star_U./nu_c;
    % WALL FUNCTION:
    if wallfunction == 1;
        u(2) = utau/Von_Karman*log(yplus(2)) + 5;
    end
    for k = 2+(wallfunction == 1):Nz+1; %if wall function used, start at 3d cell.
        %% Turbulence
        %nu_t(k) = l(k)^2*(0.5*abs((u(k+1)-u(k)))/dzc(k)+0.5*abs((u(k)-u(k-1)))/dzc(k-1));
        nu_t(k) = l(k)^2*abs((u(k+1)-u(k-1))/(dzc(k)+dzc(k-1)));
        nu(k) = nu_c+(turbulent==1)*nu_t(k);
        
        %% New a
        nu_U = (nu(k+1)*dz(k+1)+nu(k)*dz(k))/(dz(k+1)+dz(k));
        nu_D = (nu(k)*dz(k)+nu(k-1)*dz(k-1))/(dz(k)+dz(k-1));
        au = nu_U/dzc(k);
        ad = nu_D/dzc(k-1);
        ap = au+ad;
        u(k) = (u(k-1)*ad+u(k+1)*au-1/rho*dpdx(k)*dz(k))/ap;          
        
    end
    %% Enforce Boundary Conditions
    if bcswitch == 0 %velocity at both walls specified
        u(1)=-u(2)+2*uwall1;
        u(end)=-u(end-1)+2*uwall2;
    elseif bcswitch == 1 %gradient at upper boundary specified
        u(1)=-u(2)+2*uwall1;
        u(end)=u(end-1);
    elseif bcswitch == 3 %gradient at lower boundary specified
        u(1)=u(2);
        u(end)=-u(end-1)+2*uwall2;
    elseif bcswitch == 2 %velocity at upper wall, tauw at lower wall
        u(end)=-u(end-1)+2*uwall2;
        % WALL FUNCTION:
        % u(2) = utau/Von_Karman*log(yplus(2)) + 5;
    end
    Qnew = u(2:end-1)*dz(2:end-1)';
    residue = (abs(Qnew-Q))/Q*(prescribeswitch == 1)+abs(u_change-mean(u))/u_change*(prescribeswitch ==0);
    u_change = mean(u);
    iter = iter+1
end
Reynolds = u_change*H/nu_c
