clear all
close all
clc
%% INPUT VARIABLES
%% PROJECT
project        = 1;  %project = on
delta_wall     = 50; % wall boundary layer size (20 - 100 m)
u_end_boundary = 10; % wind speed at height of delta_wall (-20 - 20 m/s);
%% Steady State or transient
Steady_State_on = 1;    % Use 1 for steady state
                        % Use 0 for transient mode
Steady_State_Start = 1; % When running transient, use this to calculate the
                        % steady state first and then run transient

unsteady_function = 0;  % Use only for transient, use 0 for prescribed pressure
                        % Use 1 for flow rate
                        % These can be modified in the function files
Time_steps = 100;        % Use for transient mode
Delta_t = 10;         % Seconds between time steps
omega_unsteady = 0.01;   % Frequency (in Hz) at which the unsteady_functions change
tf = 1;                 % 0 for explicit time scheme

                        % 1 for implicit time scheme

% Mesh
Nz              = 100;  % number of cells in z-direction (should be even!)
Nx              = 1; 

H               = 500; % height of channel (100-1000 m)

L               = 1;    % length of channel

Mesh_type       = 2;    % Type of Mesh, 1 is refinement at both boundaries, 2 is refinement at only bottom
expansion_factor= 1.1;  %mesh expansion factor



% Boundary Conditions
uwall1          = 0;    % velocity at lower wall

uwall2          = 0;   % velocity at upper wall
tauw            = 0.1;  % wall shear stress
bcswitch        = 1;    % 0 if velocity is specified, 
                        % 1 if gradient at upper boundary is specified,
                        % 2 if wall shear stress is specified > USE THIS
                        % FOR WALL FUNCTION
                        % 3 if gradient at lower boundary is specified,
wallfunction   = 1;     %1 if wall function must be used (for wall shear stress)
ks             = 3;     %roughness: average height of obstacles

% Global Boundary Conditions
prescribeswitch = 0;    % 0 if pressure gradient prescribed, 
                        % 1 if flow rate prescrpibed

                        
%% Simulation

max_iter = 20000;       % Maximum iterations
min_residue = 1e-7;     % Value at which a solution is considered converged

                        % A good value is 1e-6, if this is lowered, weird
                        % results are considered converged.
% Turbulence on/off
turbulent       = 1;    % 0 if not tubulent;
                        % 1 if turbulent flow
Van_Driest_Damping_on = 1;  % 0 for no van Driest damping
                            % 1 for van Driest damping
                        
dudzwall        = 100;  % velocity gradient at the wall

wall            = 1;    % 1 if lower wall velocity gradient specified, 2 for upper (doesn't work yet)

rho             = 1.225;    % density of air
mu              = 15*10^-6;   %viscosity
nu              = mu/rho;   %viscosity

dpdx            = -0.0001;   % prescribed pressure gradient

Q               =  2;    % prescribed flow rate per area in m^2/s (2-dimensional)


show_analytic = 0;      %Use 1 to show analytic value

Mesh
Initialiser
if Steady_State_on==1
    Solver
elseif Steady_State_on==0 && Steady_State_Start==1
    Solver
    Transient_Solver
    tic
else
    Transient_Solver
    tic
end
Plotter
Running_time = Running_time + toc