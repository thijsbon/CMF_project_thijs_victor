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
uwall1          = 0;    % velocity at lower wall
uwall2          = 0;    % velocity at upper wall
bcswitch        = 0;    % 0 if velocity is specified, 
                        % 1 if gradient at upper boundary is specified,
                        % 3 if gradient at lower boundary is specified,
                        % 2 if wall shear stress is specified
% Global Boundary Conditions
prescribeswitch = 0;    % 0 if pressure gradient prescribed, 
                        % 1 if flow rate prescrpibed

                        
%% Simulation
max_iter = 20000;       % Maximum iterations
min_residue = 1e-8;     % Value at which a solution is considered converged
                        % A good value is 1e-6, if this is lowered, weird
                        % results are considered converged.
% Turbulence on/off
turbulent       = 1;    % 0 if not tubulent;
                        % 1 if turbulent flow
                        
dudzwall        = 100;  % velocity gradient at the wall
tauw            = 0.000001;  % wall shear stress
wall            = 1;    % 1 if lower wall velocity gradient specified, 2 for upper (doesn't work yet)
rho             = 1;    % density
mu              = 15*10^-6;   %viscosity
nu              = mu/rho;   %viscosity
dpdx            = -0.001;   % prescribed pressure gradient
Q               =  20;    % prescribed flow rate per area in m^2/s (2-dimensional)

show_analytic = 1;      %Use 1 to show analytic value

Mesh
Initialiser
Solver
Plotter
toc