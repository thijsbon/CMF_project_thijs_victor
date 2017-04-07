% Initialise values
%% Constants
mu = mu*ones(1,Nz+2);
nu_c = nu;              %Material property
nu = nu*ones(1,Nz+2);     %Total Viscosity
u = zeros(1,Nz+2);
dpdx = dpdx*ones(1,Nz+2);
Von_Karman = 0.41;      %Von Karman Constant
Karman_0 = 0.09;%0.085;
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