% -------------------------------------------------------------------------
% function to generate a rotation matrix that transforms coordinates from
% the fly's body frame to the wing frame given 3 euler angles
% 
% based on Tsevi/Attila's spring code
%
% INPUTS:
%   -phi, theta, eta : three euler angles (IN RADIANS) describing the
%   wing's orientation, corresponding to yaw, pitch and roll, respectively
%   (in the case of the wing, we might also call these stroke, deviation,
%   and wing pitch)
%
% OUTPUTS:
%   -rotM : rotation matrix to transform from body to wing coordinates.
%   Note: in wing coordinates, x direction lies along span, and y should
%   align along chord
% -------------------------------------------------------------------------
function rotM = body2WingRotMat(phi, theta, eta)

% write out sines and cosines to make things more compact AND switch sign
% on theta (to match Goldstein)
cph=cos(phi)    ; sph=sin(phi)   ;
cth=cos(-theta) ; sth=sin(-theta);
cet=cos(eta)    ; set=sin(eta)   ;

% compile rotation matrix
rotM = [cth*cph             ,     cth*sph          ,  -1*sth ; ...
        set*sth*cph-cet*sph , set*sth*sph+cet*cph  ,  cth*set ; ...
        cet*sth*cph+set*sph , cet*sth*sph-set*cph  , cth*cet ] ;
end