%% INPUT VARIABLES
% this is a test for github %
% now another test %
clear all; close all;
Nz              = 100;  % number of cells in z-direction inside domain (excluding ghost cells)
                        % (should be even!)
Nx              = 1;    % never becomes larger than 1?
H               = 1;    % height of channel
L               = 1;    % length of channel
uwall1          = 0;    % velocity at wall1
uwall2          = 0;    % velocity at wall2
bcswitch        = 0;    % 0 if velocity is specified, 
                        % 1 if gradient at upper boundary is specified, 
                        % 2 if wall shear stress is specified
prescribeswitch = 1;    % 0 if pressure gradient prescribed, 
                        % 1 if flow rate prescrpibed
turbulent       = 0;    % 0 if not tubulent;
                        % 1 if turbulent flow
dudzwall        = 100;  % velocity gradient at the wall
tauw            = 0.000001;  % wall shear stress
wall            = 1;    % 1 if lower wall velocity gradient specified, 2 for upper (doesn't work yet)
rho             = 1;    % density
mu              = 10^-4*ones(1,Nz); %viscosity
dpdx            = -1;   % prescribed pressure gradient
Q               =  1000;    % prescribed flow rate per area in m^2/s (2-dimensional)
exp             = 1.1;  %mesh expansion factor
%% MESH:
nonuniformmesh
if turbulent == 0; 
    %% i. for prescribed velocity at the wall:
    if prescribeswitch == 0; % pressure gradient prescribed
        if bcswitch == 0; %wall velocities specified
           pressvel 
        end
        %ii. for prescribed velocity gradient at the wall
        if bcswitch == 1; %wall gradient specified
           pressvgrad
        end
        if bcswitch == 2; %i.e. wall shear stress specified
           presstauw
       end 
    end %end pressure gradient prescribed

    %% for prescribed flow rate
    if prescribeswitch == 1;
       if bcswitch == 0; %i.e. wall velocities are specified
          flowvel
       end
       if bcswitch == 1; %i.e. gradient specified
          flowvgrad
       end   
       if bcswitch == 2; %i.e. wall shear stress specified
          flowtauw  
       end 
    end %end Q specified
end %end of non turbulent part

%% now for turbulent flow:
if turbulent == 1;
% calculate prandtl mixing length for each cell: (only works for open
% channel now)
prmixl;
 if prescribeswitch == 1;
%        if bcswitch == 0; %i.e. wall velocities are specified
%           flowvel
%        end
       if bcswitch == 1; %i.e. gradient specified
          tflowvgrad
       end   
%        if bcswitch == 2; %i.e. wall shear stress specified
%           flowtauw  
%        end 
  end %end Q specified
end




















                               