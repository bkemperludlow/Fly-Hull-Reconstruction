% -------------------------------------------------------------------
% function to get span and chord hat vectos based on euler angles
% note: wing angles are assumed to be in RADIANS
% -------------------------------------------------------------------
function [spanHat, chordHat] = getWingVecsFromAng(wingAngleMat) 
% % --------------------------------------------
% % sort out wing sign
% if contains(wingSide,'R','IgnoreCase',true)
%     wingSign = -1 ;
% elseif contains(wingSide,'L','IgnoreCase',true)
%     wingSign = +1 ;
% else
%     fprintf('Invalid wingSide entry: %s \n', wingSide)
%     keyboard
% end
% -----------------------------------
% read out angles from matrix
phi = wingAngleMat(:,1) ;
theta = wingAngleMat(:,2) ;
eta = wingAngleMat(:,3) ;

% -----------------------------------------------------
% get span hat vector
spanHat = [cos(phi).*cos(theta),sin(phi).*cos(theta), ...
    sin(theta) ] ;

% ------------------------------------------------------------
% get chord hat vector
% (based on the idea that, in the wing-tied frame of reference, the chord
% vector lies in the y'' direction. then just apply rotation matrix from
% wing to body frame)
theta_mod = -1*theta ;
%phi_mod = -1*phi ; %-1*wingSign*phi ; %

chordHat = [sin(eta).*sin(theta_mod).*cos(phi) - cos(eta).*sin(phi), ...
    sin(eta).*sin(theta_mod).*sin(phi) + cos(eta).*cos(phi) , ...
    cos(theta_mod).*sin(eta)] ;

% chordHat = [sin(eta).*sin(theta_mod).*cos(phi) + wingSign*cos(eta).*sin(phi), ...
%     sin(eta).*sin(theta_mod).*sin(phi) - wingSign*cos(eta).*cos(phi) , ...
%     cos(theta_mod).*sin(eta)] ;

% body2wing = [cth*cph            cth*sph           (-sth) ; ...
%     (sps*sth*cph-cps*sph) (sps*sth*sph+cps*cph) cth*sps ; ...
%     (cps*sth*cph+sps*sph) (cps*sth*sph-sps*cph) cth*cps ] ;
end