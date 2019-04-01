%--------------------------------------------------------------------------
% function to be used alongside largePertCorrection.m to calculate the wing
% angles, since we already have relevant vectors transformed into the
% body frame
%--------------------------------------------------------------------------
function anglesBodyFrame = calcWingAnglesLargePert(rotSpanR, rotChordR, ...
    rotSpanL, rotChordL)
%------------------------
% params 
rad2deg = (180/pi) ; 

%----------------------------------------
% calculate angles from wing vectors
phiRdeg   = atan2(rotSpanR(2),rotSpanR(1)) * rad2deg;
thetaRdeg = asin(rotSpanR(3)) * rad2deg ;
etaRdeg = calcEta(rotSpanR, rotChordR,'right') ;

phiLdeg   = atan2(rotSpanL(2),rotSpanL(1)) * rad2deg;
thetaLdeg = asin(rotSpanL(3)) * rad2deg ;
etaLdeg = calcEta(rotSpanL, rotChordL,'left') ;

%------------------------------
% store angles in array
anglesBodyFrame = [0 0 phiRdeg thetaRdeg etaRdeg phiLdeg thetaLdeg etaLdeg] ;

end