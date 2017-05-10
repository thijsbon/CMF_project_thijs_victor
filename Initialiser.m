% Initialise values
tic
Running_time = 0;
%% Constants
mu = mu*ones(1,Nz+2);
nu_c = nu;              %Material property
if Steady_State_on==1
    nu = nu*ones(Nz+2,1);     %Total Viscosity
    u = zeros(Nz+2,1);
    dpdx = dpdx*ones(Nz+2,1);
elseif Steady_State_on==0
    nu = nu*ones(Nz+2,Time_steps+1);     %Total Viscosity
    u = zeros(Nz+2,Time_steps+1);
    dpdx = dpdx*ones(Nz+2,Time_steps+1);
end
Von_Karman = 0.41;      %Von Karman Constant
Karman_0 = 0.09;%0.085;
Cs = 0.17;              %Lilly-Smagorinsky constant
Van_Driest_A = 25;      %van Driest constant
Karman_ratio = Karman_0/Von_Karman;
Boundary_Layer_Size = H-0.5*H*(bcswitch==0);
l1 = Von_Karman*zc.*(zc<Boundary_Layer_Size*Karman_ratio)+Karman_0*Boundary_Layer_Size*(zc>Boundary_Layer_Size*Karman_ratio);
if Mesh_type == 1; % symmetric l for symmetric mesh

    l = [l1(1:Nz/2+1) flip(l1(1:Nz/2+1))];
elseif Mesh_type == 2; % asymmetric l for asymmetric mesh
    l = l1;

end;

%% calculate yplus etc.
if bcswitch == 2; %i.e. tau_wall specified at bottom
   utau = sqrt(tauw/rho);
   yplus = zc.*utau./nu_c;
   if yplus(2) > 5;
       disp('y_1^+ > 5, wall function is used')
   end
end

%% DUMMY VARIABLES
dt = 1;
if unsteady_function==1
        Q_original = Q;
end
%% Solution method
if prescribeswitch == 0
    
elseif prescribeswitch == 1
    if bcswitch == 0
        
    elseif bcswitch == 1
        u = Q/sum(dz)*ones(1,Nz+2);
        dpdx = -12*Q*mu(1).*ones(1,Nz+2)/H^3;
    elseif bcswitch == 2
        
    end
end

%format shortEng
%format compact
format shortG
Time_it_takes = 0;
Running_time = toc;
