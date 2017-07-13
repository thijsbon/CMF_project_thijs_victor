
function p_at_t = unsteady_dpdx(t,Delta_t,dp_0,omega_unsteady)
%UNSTEADY_DPDX Summary of this function goes here
%   Detailed explanation goes here
    p_at_t = 1*dp_0*cos(omega_unsteady*Delta_t*t+pi/4);%-dp_0;


end

