%% FLOWVGRAD
% for prescribed flow rate and wall velocity gradient;
Qdif = 1;
       u_init = Q/sum(dz)*ones(1,Nz); %uniform velocity as initial guess;
       u = u_init; 
       dpdx = -12.*Q.*mu(1).*ones(1,Nz)/H^3; %initial dpdx, not specified at boundaries!
       count = 0;
        while max(abs(Qdif)) > 0.000001;
            u(1) = uwall1;
                for k=2:Nz-1;
                    
                    dpdx(k) = dpdx(k)+Qdif*0.01;
                    %u(k) = 0.5*(u(k-1)+u(k+1))-0.5./mu(k)*dpdx.*dz(k).^2;
                    mudown = (mu(k)*dz(k)+mu(k-1)*dz(k-1))./(dz(k)+dz(k-1));
                    muup   = (mu(k)*dz(k)+mu(k+1)*dz(k+1))./(dz(k)+dz(k+1)); %interpolate viscosity
                    u(k) = (muup/dzc(k)+mudown/dzc(k-1))^-1*(muup/dzc(k)*u(k+1)+mudown/dzc(k-1)*u(k-1)-1.*dpdx(k).*dz(k));
                    
                end
            dpdx(end) = dpdx(end)+Qdif*0.01;
            u(end) = u(end-1)+dudzwall*dzc(end)-1/mu(end)*dpdx(end)*dzc(end)*dz(end);
            Qnew = sum(u.*dz);
            Qdif = Qnew - Q;
            count = count+1;
            Qtimerate(count) = Qdif;
        end
        figure
        Qan = 0.5*dudzwall*H^2-1/(3*mu(1))*dpdx(end)*H^3;
        dpdxan = (3*mu(1)*(dudzwall/(2*H)-Q/H^3)); %analytical pressure gradient
        uan = 0.5./mu(1).*dpdxan.*(zc.^2)+zc*(dudzwall-1./mu(1)*dpdxan*H);
        %uan = zc*(dudzwall-1./mu(1)*dpdx(3)*H)
        plot(u,zc,'or',uan,zc,'b'); 