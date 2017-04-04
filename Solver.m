iter = 0;
residue = 10*min_residue;
Qnew = 0;
while iter<max_iter && residue>min_residue
    %% velocity change (for residue determination)
    u_change = mean(u);
    
    %% dpdx change (in case of the prescribed flow rate)
    dpdx = dpdx+(max_iter-iter)*mean(mu)*(Qnew-Q)/H^3*(prescribeswitch == 1);
    
    %% 
    y_star_1 = zf(2:end).*u_star./nu_c;
    y_star_2 = (H-zf(2:end)* (bcswitch==0) ).*u_star./nu_c;
    
    for k =2:Nz-1
        %% Turbulence
        % Mixing Length
        l = (min(y_star_1(k),y_star_2(k))>=35 && min(y_star_1(k),y_star_2(k))<400)*Von_Karman*min(zf,H-zf)+(min(y_star_1(k),y_star_2(k))>400)*0.085*Boundary_Layer_Size;
        u_star(k) = sqrt(abs((0.5*(u(k+1)-u(k))/dz(k)+0.5*(u(k)-u(k-1))/dz(k-1))/rho));
        nu_t(k) = l(k)^2*u_star(k)^2;
        nu(k) = nu_c+(turbulent==1)*nu_t(k);
        
        %% New a
        au = (nu(k+1)+nu(k))/dzc(k);
        ad = (nu(k)+nu(k-1))/dzc(k-1);
        ap = au+ad;        
        u(k) = (u(k-1)*ad+u(k+1)*au-1/rho*dpdx(k)*dz(k))/ap;          
        
    end
    %% Enforce Boundary Conditions
    if bcswitch == 0
        u(1)=-u(2)+2*uwall1;
        u(Nz)=-u(Nz-1)+2*uwall2;
    elseif bcswitch == 1
        u(1)=-u(2)+2*uwall1;
        u(Nz)=u(Nz-1);
    elseif bcswitch == 3
        u(1)=u(2);
        u(Nz)=-u(Nz-1)+2*uwall2;
    end
    Qnew = u*dz';
    residue = (abs(Qnew-Q))/Q*(prescribeswitch == 1)+abs(u_change-mean(u))/u_change;
    iter = iter+1
end
Qnew
