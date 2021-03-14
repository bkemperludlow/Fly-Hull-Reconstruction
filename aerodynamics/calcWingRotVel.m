% -------------------------------------------------------------------------
% function to calculate wing rotational velocity, as measured in both wing
% and body frames (in body frame, we calculate wing rotation w.r.t. body)
%
% phi, theta, and psi are the wing stroke, deviation, and rotation angles,
% respectively. *_dot is that quantity's derivative w.r.t. time
% -------------------------------------------------------------------------
function [omegaWingBody, omegaWingWing] = calcWingRotVel(phi, theta, psi, ...
    phi_dot, theta_dot, psi_dot)
% ----------------------------------------------------------------------
%% write out sines and cosines 
% make things more compact AND switches sign on theta (to match Goldstein)
cph=cos(phi)    ; sph=sin(phi)   ;
cth=cos(-theta) ; sth=sin(-theta);
cps=cos(psi)    ; sps=sin(psi)   ;

modThetadot    = -theta_dot ;

% ----------------------------------------------------------------------
%% wing angular velocity as measured in the wing frame of reference.
omegaWingWing = [ psi_dot - phi_dot.*sth   , ...
    modThetadot.*cps  + phi_dot.*cth.*sps  , ...
    -modThetadot.*sps + phi_dot.*cth.*cps ] ;

% ----------------------------------------------------------------------
%% wing angluar velocity as measured in the fly's body frame
omegaWingBody = [psi_dot.*cth.*cph - modThetadot.*sph , ...
    psi_dot.*cth.*sph  + modThetadot.*cph , ...
    phi_dot - psi_dot.*sth] ;

end