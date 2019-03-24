function [f fLift fDrag fRot] = quasiSteadtForceRotation(vbs, geom)
% quasi steadt aerodynamic force including a translational U^2 force and
% a rotational force, as formulated by Sane and Dickinson, JEB 2002
% (Revised quasi-steady model)

% vbs is the wing tip velocity in the wing frame of reference
% angular velocities should be given in rad/sec

% the "translational" force is the same as in quasiSteadyForceSimple.m

rv = zeros(4, 1) ;
rv(1) = geom.span / 2;
rv(4) = 1 ;

vmat = wedge(vbs) * rv ;
%vxp = vmat(1,:) ;
vyp = vmat(2) ;
vzp = vmat(3) ;

% absv = sqrt(vyp.^2 + vzp.^2) ;
U     = [ 0 vyp vzp] ;

U2    = (vyp.^2 + vzp.^2) ;
if (U2==0)
    %disp('zero velocity deal with this here. return zero force') ;
    f = zeros(3,1) ;
    fLift = zeros(3,1) ;
    fDrag = zeros(3,1) ;
    fRot  = zeros(3,1) ;
    return ;
end

Uhat  = U / sqrt(U2) ;
alpha = atan2(vzp, vyp) ; 

airDensity = geom.airDensity ; 

r22 = geom.r22 ; % 0.4 ; % 0.3 ; % geometric factor, second moment of area
wingArea = geom.span * geom.chord * pi / 4; % area of an ellipse

CLMAX = 1.8 ; % see Dickinson Science 1999 and Whitney & Wood JFM 2010
CDMAX = 3.4 ;
CD0   = 0.4 ;
CROT  = 1.55 ; 
rotShapeIntegral = 0.5 ; % estimation

% lift and drag coefficients
%CD = (CDMAX/2) * ( 1 - cosd(2*alpha)) ;
CD = (CDMAX+CD0)/2 - ((CDMAX-CD0)/2) * cos(2*alpha) ;
CL = CLMAX * abs(sin(2*alpha)) ;

fDrag = -0.5 * CD * airDensity * wingArea * r22 * Uhat * U2 ;


spanHat = [1 0 0 ] ; 
Lhat = cross(Uhat, spanHat) ; % perpendicular to Uhat
Lhat = Lhat / norm(Lhat) ; % but still need to determine sign

% Lhat and Uhat should be on different sides on the wing, so they should
% have different signs of their zprime component

% sign(Lhat(3)) * sign(Uhat(3)) is negative if Lhat is on the right
% direction. Hence we have:

Lsign = -sign(Lhat(3)) * sign(Uhat(3)) ;

Lhat  = Lhat * Lsign ;

fLift = + 0.5 * CL * airDensity * wingArea * r22 * Lhat * U2 ;

fDrag = fDrag' ;
fLift = fLift' ;

% ROTATIOAL FORCE
% ---------------

psidot = vbs(4) ; % x-axis in the wing frame is the span direction

fRotMag = CROT * airDensity * sqrt(U2) * psidot * ...
        geom.chord^2 * geom.span * rotShapeIntegral ;

% if psidot>0 then force is directed along zprime axis
fRot = [0 ; 0 ; fRotMag ] ;

f = fDrag + fLift + fRot;

    
%{
if (includeRotationalForce)
    
    % rotational force
    % ----------------
    % see eq. 11 in Sane and Dickinson JEB 2001
    % calculate magnitude and direction separately
    
    currt = t(it) ;
    etadot_rad_sec = eval_fd(currt, fdEtaR_lab, 1) * (pi/180) * 1000 ; % 1000 from msec to sec
    % NEED TO USE OMEGA=psidot + phidot*sin(theta)
    
    psiDot = 
    
    spanHat = data.rightSpanHats(it,:);
    rotHat = cross(spanHat, chordHat) * sign(etadot_rad_sec) ;
    rotHat = rotHat/norm(rotHat) ;
    
    frot = CROT * airDensity * sqrt(U2) * abs(etadot_rad_sec) * ...
        wingChordMeters^2 * wingSpanMeters * rotShapeIntegral ;
    
    FrotR(it,:) = frot * rotHat ;
    
    
end

%}
return



