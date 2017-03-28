%% PRESSVGRAD
% for prescribed pressure gradient and wall velocity gradient
u_old = ones(1,Nz);
        u = [0 zeros(1,Nz-2) 0];
        test = 1;
        count = 0;
        while max(test) > 0.000001;
            u(end) = u(end-1) + dudzwall*dzc(end)-1/mu(end)*dpdx*dzc(end)*dz(end);
            u(1) = uwall1;
                for k=2:Nz-1;
                    %u(k) = 0.5*(u(k-1)+u(k+1))-0.5./mu(k)*dpdx.*dz(k).^2;
                    mudown = (mu(k)*dz(k)+mu(k-1)*dz(k-1))./(dz(k)+dz(k-1));
                    muup   = (mu(k)*dz(k)+mu(k+1)*dz(k+1))./(dz(k)+dz(k+1)); %interpolate viscosity
                    % in the Kenjeres way:
                    ad = mudown/dzc(k-1);
                    au = muup/dzc(k);
                    ak = ad + au;
                    b = -dpdx*dz(k);
                    u(k) = 1/ak*(ad*u(k-1)+au*u(k+1)+b);
                    %u(k) = (muup/dzc(k)+mudown/dzc(k-1))^-1*(muup/dzc(k)*u(k+1)+mudown/dzc(k-1)*u(k-1)-1.*dpdx.*dz(k));
                    %u(k) = (mu(k)/dzc(k)+mu(k-1)/dzc(k-1))^-1*(mu(k)/dzc(k)*u(k+1)+mu(k-1)/dzc(k-1)*u(k-1)-1.*dpdx.*dz(k));
                end
            test = abs(u-u_old);
            u_old = u;
            count = count+1;
        end
        figure
        Q = sum(u.*dz);
        Qan = 0.5*dudzwall*H^2-1/(3*mu(1))*dpdx*H^3;
        dpdxan = (3*mu(1)*(dudzwall/(2*H)-Q/H^3));
        uan = 0.5./mu(1).*dpdx.*(zc.^2)+zc*(dudzwall-1./mu(1)*dpdx*H);
        
        plot(uan,zc,'-b',u,zc,'or');