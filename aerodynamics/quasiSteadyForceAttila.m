function f = quasiSteadyForceAttila(vbs, geom)

% constants
CT    = 1.833;
CR    = 1.5  ;
CD0   = 0.21 ;
CDPI2 = 3.35 ;
MU1   = 0.2  ;
MU2   = 0.2  ;

RHOF  = geom.airDensity ;% air density in kg/m^3

Nblades  = 10   ;

rx = geom.span / 2 ;
ry = geom.chord / 2 ;

r = linspace(-1, 1, Nblades) * rx ;
rv = zeros(4, Nblades) ;
rv(1,:) = r ;
rv(4,:) = 1 ;

vmat = wedge(vbs) * rv ;

%vxp = vmat(1,:) ;
vyp = vmat(2,:) ;
vzp = vmat(3,:) ;

absv = sqrt(vyp.^2 + vzp.^2) ;
alpha = atan2(vzp, vyp) ;

% local chord of an ellipse
chord = 2 * ry * sqrt( 1 - (r/rx).^2) ;

gamma = -CT * chord .* sin(2*alpha) .* absv / 2 + ...
    CR * vbs(4) * chord.^2 / 2 ;

fnu = RHOF * chord / 2 .* ( CD0*cos(alpha).^2 + CDPI2*sin(alpha).^2 ) .* absv ;

taunu = pi * RHOF * (chord/2).^4 * (MU1*geom.freq + MU2*abs(vbs(4))*vbs(4)) ;

circy = -RHOF * gamma .* vzp ;
circz =  RHOF * gamma .* vyp ;

dragy = -fnu .* vyp ;
dragz = -fnu .* vzp ;

fyp = circy + dragy ;
fzp = circz + dragz ;

% combine to form the wrench along a blade element
wb = [ zeros(1,Nblades) ; fyp ; fzp ; -taunu ; -r.*fzp ; r.*fyp] ;

% integrate
dr = r(2) - r(1) ;
f = trapz(wb,2) * dr ;
    

%{
        # wrench along a blade element
        wb = _n.array([_n.zeros(fyp.shape), fyp, fzp, -taunu, -r*fzp, r*fyp])
        fbs[0, :] = _n.trapz(wb, r)
        
        return fbs
%}

return

