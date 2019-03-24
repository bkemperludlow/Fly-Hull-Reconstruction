function [wrenchWingFrame forceLabFrame wing2body body2wing torqueLabFrame vlab] = ...
    quasiSteadyForce(phi, theta, psi, ...
    phidot, thetadot, psidot, vBody,  vCent, geom, modelType)

% phi, theta, psi: wing Euler angles 
% phidot, thetadot, psidot:  angular velocities in the LAB frame of reference
% vBody is the velocity of the fly's center-of-mass (3 vector)
%
% IMPORTANT: the angle theta and thetadot inputs are defined such that a
% positive theta takes the wing ABOVE the xy plane. This is the opposite
% direction of the definition of Goldstein yaw-pitch-roll, that take theta
% as a right handed rotation along the new y-axis (after phi rotation about
% the lab z-axis). So, anywhere we're using here the Goldstein notaion we
% first transform theta --> (-theta)
%
% vCent (OPTIONAL) is the velocity of the wing center (3 vector) 
% it is OPTIONAL - if vCent is not specified, then it is calculated from
% the angular velocity of the wing (w x rcm) to which we add vBody
%
% geom is a structure with the following fields:
%       geom.span  = span ;
%       geom.chord = chord ;
%       geom.thick = chord / 50 ;
%       geom.mass  = 1e-6 * 0.003 ; % wing mass ;
%       geom.freq  = wingFreq ;
% modelType is;
%    2: for simple quasi steady with U^2 "translational" force only
%    1: "Attila's" quasi steady model
%    3: simple QS with rotational force
%}

cph=cos(phi)    ; sph=sin(phi)   ;
cth=cos(-theta) ; sth=sin(-theta);
cps=cos(psi)    ; sps=sin(psi)   ;

body2wing = [cth*cph            cth*sph           (-sth) ; ...
    (sps*sth*cph-cps*sph) (sps*sth*sph+cps*cph) cth*sps ; ...
    (cps*sth*cph+sps*sph) (cps*sth*sph-sps*cph) cth*cps ] ;

wing2body = body2wing' ;

% Goldstein appendix = attila's thesis with -1*theta
%  should have used the name omegaWing rather than omegaBody
modThetadot=-thetadot ;
omegaWing = [ psidot - phidot*sth   ; ... % note minus sign is definition of sth above, so minus here is ok
    modThetadot*cps  + phidot*cth*sps  ; ...
    -modThetadot*sps + phidot*cth*cps ] ;

% 6-vector velocity in the body frame of reference (bs=body system)

omegaLab = [psidot*cth*cph - modThetadot*sph ; ...
    psidot*cth*sph  + modThetadot*cph ; ...
    phidot - psidot*sth] ; % XXX


if (isempty(vBody))
    vBody = 0 ;
end

if (isempty(vCent))
    % http://planning.cs.uiuc.edu/node102.html
    % rotation matrices
    %{
    R1 = [ cph -sph 0 ; ...
        sph cph 0 ; ...
        0    0  1 ] ;
    
    R2 = [cth 0 sth ; ...
        0     1  0 ; ...
        -sth  0  cth] ;    %
    %}
    rcent = (wing2body)*[geom.span/2 ;-geom.chord/2 ;0] ;     % wing center position
    vlab = cross(omegaLab, rcent) + vBody; % velocity of wing center in the lab frame
else
    vlab = vCent ;
end

vWingFrame = body2wing*vlab ; % wing velocity in the wing frame (vBody is included here)

vbs = [vWingFrame ; omegaWing ];

switch (modelType)
    case 1
        wrenchWingFrame = quasiSteadyForceAttila(vbs, geom) ;
    case 2
        [f,~,~] = quasiSteadtForceSimple(vbs, geom) ; 
        wrenchWingFrame = [f ;0 ;0 ;0] ;
    case 3
        [f,~,~,~] = quasiSteadtForceRotation(vbs, geom) ; 
        wrenchWingFrame = [f ;0 ;0 ;0] ;
end

forceLabFrame  = wing2body * wrenchWingFrame(1:3) ;
torqueLabFrame = wing2body * wrenchWingFrame(4:6) ;

% NO. output an additional parameter from Attila's force - the wb matrix

% wrenchWingFrame(4:6) from attila's code gives the torque on the wing
% center expressed in the wing frame of reference.
% the question is how to transform it to the body CM in the lab frame of
% reference. See eq. 2.66 in M.L.S. Book. Do we need to consider the bottom
% left term in the matrix? Has it been considered when taking r.*fyp etc...
% ??? YES WE NEED TO CONSIDER ALL TERM


return
