%% PRMIXL
% calculate prandtl mixing length stuff

% !!! only works for open channel now !!!

kappa = 0.4; %Von Karmann Constant
Ctr   = 0.2; %transition at ztr = Ctr*delta
lambda = 0.2; % above ztr, l = lambda*delta
%calculate height of boundary layer delta:
if prescribeswitch == 0;
   delta = 0.5*H;
  
end
if prescribeswitch == 1;
   delta = H;
    %calculate prandtl mixing length prl:
    ztr = Ctr*delta; %calculate transition z;
    prl = zeros(length(zc),1); % !!!!!!!!!!! prl defined at cell centers, like mu. 
    prl(zc<ztr) = kappa*zc(zc<ztr);
    prl(zc>ztr) = lambda*zc(zc>ztr);
end

