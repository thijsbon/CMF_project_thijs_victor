% Initialise values
%% Constants
mu = mu*ones(1,Nz);
u = zeros(1,Nz);
dpdx = dpdx*ones(1,Nz);

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