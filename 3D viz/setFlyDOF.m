function setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, ...
                   bodyRcm, bodyYPR, rightYPR , leftYPR, thetab0, dL) 

% angles in deg


wingYawAxis = [sin(thetab0) 0 cos(thetab0) ] ;
wingThetaAxis = [cos(thetab0) 0 -sin(thetab0) ] ;

Trw = makehgtform('translate',[0 -dL 0]) ;
rightWingRoll  = makehgtform('yrotate',thetab0-rightYPR(3)) ;
rightWingPitch = makehgtform('axisrotate',wingThetaAxis, -rightYPR(2) ) ;
%disp('zero and sign of right wing yaw might need to change here') ;
rightWingYaw = makehgtform('axisrotate',wingYawAxis, pi/2-rightYPR(1) ) ;
set(rightWingGrp,'Matrix',Trw*rightWingYaw*rightWingPitch*rightWingRoll) ;


Tlw = makehgtform('translate',[0 +dL 0]) ;
leftWingRoll  = makehgtform('yrotate',thetab0-leftYPR(3)) ;
leftWingPitch = makehgtform('axisrotate',wingThetaAxis, +leftYPR(2) ) ;

%disp('zero and sign of left wing yaw might need to change here') ;
leftWingYaw = makehgtform('axisrotate',wingYawAxis, -pi/2+leftYPR(1) ) ;
set(leftWingGrp,'Matrix',Tlw*leftWingYaw*leftWingPitch*leftWingRoll) ;


% IMPORTANTOS: rotate the whole fly only AFTER wing angles have been set
% test yaw, pitch, roll of whole fly
Rroll  = makehgtform('xrotate',bodyYPR(3)) ;
Ryaw   = makehgtform('zrotate',bodyYPR(1)) ;
Rpitch = makehgtform('yrotate',-bodyYPR(2)) ;
Rtrans = makehgtform('translate',bodyRcm) ;

M = eye(4) ; % get(flyGrp,'Matrix') ;
set(flyGrp,'Matrix',Rtrans*Ryaw*Rpitch*Rroll*M) ;

end