if Mesh_type == 1
    % builds 1D vector with finer mesh near walls
    dz1     = 0.5*H/(1+sum(exp.^(1:Nz/2-1))); %size of first cell
    dzvec   = [dz1 dz1 exp.^(1:Nz/2-1)*dz1];  %vector with sizes of cells in lower half of domain, ghost cell has same size as first cell in domain
    % for environmental flow with open upper wall: perhaps only grid refinement
    % near bottom?
    zf      = [0 cumsum(dzvec(2:end))];       %vector with coodinates of cell faces in lower half of domain (starts at 2 s.t. ghost cell is not taken into account) 
    zf      = [zf H-flip(zf(1:end-1))];       %extend to upper half of domain
    zf      = [-zf(2) zf H+(H-zf(end-1))];    %add faces of ghost cells (same size as first and last cell)                                       
    dz      = diff(zf);                       %size of control volumes
    zc      = zf(1:end-1)+dz/2;               %center coordinates
    dzc     = diff(zc);                       %differences between centers;

    %plot(0,zf,'xr',0,zc,'ob'); %legend('faces','centers'); %legend is messed up? 
end