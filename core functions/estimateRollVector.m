function [rhoTimes, rollVectors] = estimateRollVector(data)
%NEED TO ADJUST TO BE CALCULATING ANGLES IN PLANE PERPENDICULAR TO AHAT

%If rhoTimes and rollVectors are not defined (when no manual 
%correction has been performed) this guesses them 
tol_1 = -.7 ; %-.9
%tol_2 = .3 ;  %.2

AHat = data.AHat ;
wingTipR = data.rightWingTips ;
wingTipL = data.leftWingTips ;
wingCMR = data.rightWingCM ;
wingCML = data.leftWingCM ;
bodyCM = data.bodyCM ;
spanHatR = data.rightSpanHats ;
spanHatL = data.leftSpanHats ;

Nimages = length(AHat) ;

wingTipR_bodyFrame = wingTipR - bodyCM ;
wingTipL_bodyFrame = wingTipL - bodyCM ;

wingCMR_bodyFrame = wingCMR - bodyCM ;
wingCML_bodyFrame = wingCML - bodyCM ;


Ut_R = diff(wingTipR_bodyFrame) ;
Ut_R_xy = [Ut_R(:,1:2) zeros(length(Ut_R),1)] ;
Ut_R_xy = Ut_R_xy ./ repmat(myNorm(Ut_R_xy),1,3) ;
Ut_R_xy = [0, 0, 0; Ut_R_xy] ;

Ut_L = diff(wingTipL_bodyFrame) ;
Ut_L_xy = [Ut_L(:,1:2) zeros(length(Ut_L),1)] ;
Ut_L_xy = Ut_L_xy ./ repmat(myNorm(Ut_L_xy),1,3) ;
Ut_L_xy = [0, 0, 0; Ut_L_xy] ;

AHat_xy = [AHat(:,1) AHat(:,2) zeros(length(AHat),1)] ;
AHat_xy = AHat_xy ./ repmat(myNorm(AHat_xy),1,3) ;

spanHatR_perpA = spanHatR - repmat(dot(spanHatR,AHat,2),1,3).*AHat ;
spanHatR_perpA = spanHatR_perpA ./ repmat(myNorm(spanHatR_perpA),1,3) ;
spanHatL_perpA = spanHatL - repmat(dot(spanHatL,AHat,2),1,3).*AHat ;
spanHatL_perpA = spanHatL_perpA ./ repmat(myNorm(spanHatL_perpA),1,3) ;

%goodRollPoints are frames during which 1) the projection of the body
%axis into the xy plane is roughly antiparallel to the wing tip
%velocity (also projected into the xy plane) AND 2) the projections of
%the spans into the xy plane are roughly perpendicular to the
%xy-projected body axis
goodRollPoints = find((dot(AHat_xy,Ut_L_xy,2) < tol_1) & (dot(AHat_xy,Ut_R_xy,2) < tol_1) ...
    & ( dot(spanHatR_perpA,spanHatL_perpA,2) < tol_1)) ;

% rhoTimes defined in FRAMES (i think that's right)
rollVectors = zeros(Nimages,3) ;
if ~isempty(goodRollPoints)
    rhoTimes = goodRollPoints' ;
    %rollvector goes from R --> L wing (normally hinge, in this case CoM)
    rollVec = wingCML_bodyFrame(goodRollPoints,:) - wingCMR_bodyFrame(goodRollPoints,:) ;
    %need to make sure rollHats orthogonal to AHat!
    rollVec = rollVec ./ repmat(myNorm(rollVec),1,3) ;
    
    rollVec_proj = rollVec - repmat(dot(AHat(goodRollPoints,:),rollVec,2),1,3).*AHat(goodRollPoints,:) ;
    rollVec_proj = rollVec_proj ./ repmat(myNorm(rollVec_proj),1,3) ;
    
    for j = 1:length(goodRollPoints)
        rollVectors(goodRollPoints(j),:) = rollVec_proj(j,:) ;
    end    
else
    rhoTimes = 1:Nimages ;
end


end