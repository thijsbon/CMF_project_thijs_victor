% Initialise values
%% Constants
mu = mu*ones(1,Nz);
nu_c = nu;              %Material property
nu = nu*ones(1,Nz);     %Total Viscosity
u = zeros(1,Nz);
dpdx = dpdx*ones(1,Nz);
Von_Karman = 0.41;      %Von Karman Constant
Boundary_Layer_Size = H-0.5*H*(bcswitch==0);
u_star = zeros(1,Nz);

%% DUMMY VARIABLES
dt = 1;

%% Solution method
if prescribeswitch == 0
    
elseif prescribeswitch == 1
    if bcswitch == 0
        
    elseif bcswitch == 1
        u = Q/sum(dz)*ones(1,Nz);
        dpdx = 12*Q*mu(1).*ones(1,Nz)/H^3;
    elseif bcswitch == 2
        
    end
end