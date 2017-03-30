if Mesh_type == 1
    % NONUNIFORMMESH
    % builds 1D vector with finer mesh near walls
    dz1     = 0.5*H/(0.5+sum(exp.^(1:Nz/2-1)));
    dz2     = [0.5*dz1 exp.^(1:Nz/2-1)*dz1];
    zf      = cumsum(dz2);
    zf      = [zf H-flip(zf(1:end-1))];
    zf      = [-zf(1) zf H+zf(1)]; % face coordinates
    dz      = diff(zf); % size of control volumes
    zc      = zf(1:end-1)+dz/2; % center coordinates
    dzc     = diff(zc); % differences between centers;

    %plot(0,zf,'xr',0,zc,'ob')
end