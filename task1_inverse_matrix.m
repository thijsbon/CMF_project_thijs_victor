%% i. for prescribed velocity at the wall:
% vecta = ones(length(zc)-2,1); %vector with only ones
% vectb = -2.*vecta;            %vector with only -2
% vectc = vecta;                %vector with only ones
% vecta = [0; vecta];           %make first value 0
% vectb = [1 ;vectb; 1];        %make first value one
% vectc = [vectc ;0];           %make last value 0
% A = diag(vectb)+diag(vecta,1)+diag(vectc,-1); %build matrix A
% vectd = 1./mu.*dpdx.*dz.^2;
% vectd(1) = uwall1; vectd(end) = uwall2;
%% ii. for prescribed velocity gradient at the wall:
% vecta = ones(length(zc)-2,1); %vector with only ones
% vectb = -2.*vecta;            %vector with only -2
% vectc = vecta;                %vector with only ones
% vecta = [0.5; vecta];           %make first value 0
% vectb = [-0.5 ;vectb; 0.5];        %make first value one
% vectc = [vectc ;-0.5];           %make last value 0
% A = diag(vectb)+diag(vecta,1)+diag(vectc,-1); %build matrix A
% vectd = 1./mu.*dpdx.*dz.^2;
% vectd(1) = dudzwall*dz(end); vectd(end) = dudzwall*dz(end);