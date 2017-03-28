% TASK 1. 
% Develop a routine to compute the velocity in a steady fully-developed single-phase
% channel flow with variable viscosity, using the finite-volume method. Assume that the
% viscosity profile is a prescribed input to the routine. Consider three possible wall boundary
% conditions: (i) prescribed velocity at the wall, (ii) prescribed velocity-gradient at the wall,
% and (iii) prescribed relation between the shear-stress at the wall and the velocity field.
% Consider two possible global boundary conditions: (i) prescribed pressure-gradient and
% (ii) prescribed flow-rate.
%%  Input
clear all; close all;
Nz = 100; % number of cells
Nx = 1; 
H = 1; % height of channel
L = 1; % length of channel
uwall1 = 0; %velocity at wall1
uwall2 = 0; %velocity at wall2
bcswitch = 0; % 0 if velocity is specified, 1 if gradient at wall is specified
prescribeswitch = 0; % 0 if pressure gradient prescribed, 1 if flow rate prescrpibed
dudzwall = 100; %velocity gradient at the wall
rho = 1; %density
mu = 10^-4*ones(1,Nz); %viscosity
dpdx = -1; % prescribed pressure gradient
Q = 1; % prescribed flow rate (in m^2/s) (per unit depth)
%% Mesh: (UNIFORM!!)
z = linspace(0,H,2*Nz-1); % create equidistant grid
x = linspace(0,L,2*Nx-1);
zc = z(1:2:end);          %center coordinates are odd elements
zf = z(2:2:end);          %face coordinates are even elements
zf = [-zf(1) zf zc(end)+zf(1)];   %make face coordinates outside domain
dz = diff(zf); %size of control volumes
dzc = diff(zc);
%plot(zf,0,'xb',zc,0,'or')


%% i. for prescribed velocity at the wall:
if prescribeswitch == 0;
if bcswitch == 0;
u_old = ones(1,Nz);
u = [uwall1 zeros(1,Nz-2) uwall2];
test = 1;
count = 0;
while max(test) > 0.0000001;    
    for k=2:Nz-1;
        %u(k) = 0.5*(u(k-1)+u(k+1))-0.5./mu(k)*dpdx.*dz(k).^2;
        u(k) = (1/dzc(k)+1/dzc(k-1))^-1*(1/dzc(k)*u(k+1)+1/dzc(k-1)*u(k-1)-1./mu(k).*dpdx.*dz(k));
    end
test = abs(u-u_old);
u_old = u;
count = count+1;
end
upois = 1/(2*mu(1))*dpdx*(zc.^2-H*zc);
plot(u,zc,'or',upois,zc,'b-')
end
Qtest = sum(u);
%% ii. for prescribed velocity gradient at the wall
if bcswitch == 1; 
u_old = ones(1,Nz);
u = [0 zeros(1,Nz-2) 0];
test = 1;
count = 0;
while max(test) > 0.000001;
    %u(1) = u(2)-dudzwall*dz(1)-0.5*1/mu(1)*dpdx*dz(1)^2;
    u(1) = u(2)-dudzwall*dz(1)-0.5*1/mu(1)*dpdx*dzc(1)*dz(1);
    u(end) = 0;
    for k=2:Nz-1;
        u(k) = 0.5*(u(k-1)+u(k+1))-0.5./mu(k)*dpdx.*dz(k).^2;
    end
    
test = abs(u-u_old);
u_old = u;
count = count+1;
end
uan = 0.5./mu(1).*dpdx.*(zc.^2-H^2)+dudzwall.*(zc-L);
plot(uan,zc,'-b',u,zc,'or');
end
end

%% Prescribed flow rate:
if prescribeswitch == 1;
   
end

