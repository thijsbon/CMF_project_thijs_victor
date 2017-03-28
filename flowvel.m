%% FLOWVEL
% for prescribed flow rate and wall velocity
        Qdif = 100;
        u_init = Q/sum(dz)*ones(1,Nz); %uniform velocity as initial guess;
        u_init(1) = uwall1; u_init(end) = uwall2; %set wall velocities to uwall
        dpdx = -12.*Q.*mu(1).*ones(1,Nz)/H^3; %initial guess dpdx, not specified at boundaries!
        count = 0;
        u = u_init;
        while max(abs(Qdif)) > 0.00001;
            for k = 2:Nz-1;
                %calculate new pressure gradient:
                dpdx(k) = dpdx(k)+Qdif*0.01;
                %calculate new u
                u(k) = (mu(k)/dzc(k)+mu(k-1)/dzc(k-1))^-1*(mu(k)/dzc(k)*u(k+1)+mu(k-1)/dzc(k-1)*u(k-1)-1.*dpdx(k).*dz(k));
            end;
            %calculate Qnew
            pQnew = sum(u.*dz);
            Qdif = (Qnew-Q); %calculate difference between calculated Q and specified Q;
            count = count+1;
        end
        upois = -6*Q/H^3*(zc.^2-H*zc);
        plot(u,zc,'ro',upois,zc,'b')
    