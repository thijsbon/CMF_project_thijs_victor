clear all
close all
clc
tic
%% INPUT VARIABLES
% Mesh
Nz              = 100;  % number of cells in z-direction (should be even!)
Nx              = 1; 
H               = 1;    % height of channel
L               = 1;    % length of channel
Mesh_type       = 1;    % Type of Mesh, 1 is nonuniform
exp             = 1.1;  %mesh expansion factor
% Boundary Conditions
uwall1          = 0;    % velocity at wall1
uwall2          = 0;    % velocity at wall2
bcswitch        = 1;    % 0 if velocity is specified, 
                        % 1 if gradient at upper boundary is specified, 
                        % 2 if wall shear stress is specified
% Global Boundary Conditions
prescribeswitch = 1;    % 0 if pressure gradient prescribed, 
                        % 1 if flow rate prescrpibed

                        
%% Simulation
max_iter = 20000;       % Maximum iterations
min_residue = 0.001;    % Value at which a solution is considered converged
% Turbulence on/off
turbulent       = 0;    % 0 if not tubulent;
                        % 1 if turbulent flow
                        
dudzwall        = 100;  % velocity gradient at the wall
tauw            = 0.000001;  % wall shear stress
wall            = 1;    % 1 if lower wall velocity gradient specified, 2 for upper (doesn't work yet)
rho             = 1000;    % density
mu              = 10^-6;   %viscosity
dpdx            = -1;   % prescribed pressure gradient
Q               =  1000;    % prescribed flow rate per area in m^2/s (2-dimensional)

Mesh
Initialiser
Solver

toc