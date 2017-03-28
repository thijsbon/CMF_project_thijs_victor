%% PRESSVEL
% for prescribed pressure gradient and velocity at the wall
u_old = ones(1,Nz);
        u = [uwall1 zeros(1,Nz-2) uwall2];
        test = 1;
        count = 0;
            while max(test) > 0.000001;
                for k=2:Nz-1;
                    %u(k) = 0.5*(u(k-1)+u(k+1))-0.5./mu(k)*dpdx.*dz(k).^2;
                    %u(k) = (mu(k)/dzc(k)+mu(k-1)/dzc(k-1))^-1*(mu(k)/dzc(k)*u(k+1)+mu(k-1)/dzc(k-1)*u(k-1)-1.*dpdx.*dz(k));
                    mudown = (mu(k)*dz(k)+mu(k-1)*dz(k-1))./(dz(k)+dz(k-1));
                    muup   = (mu(k)*dz(k)+mu(k+1)*dz(k+1))./(dz(k)+dz(k+1)); %interpolate viscosity
                    u(k) = (muup/dzc(k)+mudown/dzc(k-1))^-1*(muup/dzc(k)*u(k+1)+mudown/dzc(k-1)*u(k-1)-1.*dpdx.*dz(k));
                end
                    test = abs(u-u_old);
                    u_old = u;
                    count = count+1;
                end
        z = linspace(0,H,2*Nz-1);
        upois = 1./(2.*mu).*dpdx.*(zc.^2-H*zc); %analytic solution
        Qpois = -1./(12.*mu).*dpdx.*H^3; %analytic flow rate
        Qtest = sum(u.*dz); %simulation flow rate
        figure
        plot(u,zc,'or',upois,zc,'b-')