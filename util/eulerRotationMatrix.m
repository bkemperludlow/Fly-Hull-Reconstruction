function M = eulerRotationMatrix(phiRad, thetaRad, psiRad)
% M = eulerRotationMatrix(phiRad, thetaRad, psiRad)
% Euler rotation matrix, where angles are given in RADIANS
% the rotation matrix represents the rotation about the given Euler angles,
% such that the rotation about theta>0 is ABOVE the horizon. This is opposite than
% Goldstein's definitions. His definitions make sense since a torque along
% positive (intermediate) y increases theta. This definition makes sense
% because this is how we're use to think about flies and airplanes.
% Note the usage of M and its transpose:
%
% M * v = a 
% a are the coordinates of the lab-frame-vector v as seen in the rotated frame
%
% M' * u = b
% b are the lab-frame coordinates of the lab-frame-vector u RoTATED by (phi, theta, psi)

thetaRad = - thetaRad ;

cph=cos(phiRad)    ; sph=sin(phiRad)   ;
cth=cos(thetaRad)  ; sth=sin(thetaRad) ;
cps=cos(psiRad)    ; sps=sin(psiRad)   ;

M = [cth*cph                    cth*sph         (-sth) ; ...
    (sps*sth*cph-cps*sph) (sps*sth*sph+cps*cph) cth*sps ; ...
    (cps*sth*cph+sps*sph) (cps*sth*sph-sps*cph) cth*cps ] ;

return
