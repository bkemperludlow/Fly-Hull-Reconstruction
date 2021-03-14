% -------------------------------------------------------------------------
% new attempt at estimating roll vectors from data, which we can later use
% to calculate roll angle for full movie
%
% If rhoTimes and rollVectors are not defined (when no manual 
% correction has been performed) this guesses them 
% -------------------------------------------------------------------------
function [rhoTimes, rollVectors] = ...
    estimateRollVector(data, largePertFlag, interpFlag)
% ---------------------------
%% params and inputs
if ~exist('largePertFlag','var') || isempty(largePertFlag)
   try
       largePertFlag = checkAHat(data) ; 
   catch
       largePertFlag = false ; 
   end
end
if (~exist('interpFlag','var') || isempty(interpFlag)) && ~largePertFlag
   interpFlag = false ; 
elseif (~exist('interpFlag','var') || isempty(interpFlag)) && largePertFlag
    interpFlag = true ; 
end

% tolerances for dot product check
tol_1 = -.7 ; %-.9
%tol_2 = .3 ;  %.2

% --------------------------------------------
%% calculate angles to align AHat with x axis
% this should make left/right checks easier
if largePertFlag
    % in cases where gimbal lock may be a problem, use a frame-by_frame
    % estimate for the body pitch and yaw angles
    [~, ~, rotM_YP] = calcPitchLargePert(data) ;
else
    % otherwise proceed as normal
    rawPsi  = zeros(data.Nimages,1) ; % yaw
    rawBeta = zeros(data.Nimages,1) ; % pitch
    rotM_YP = zeros(3,3,data.Nimages) ; % rotation matrices that will unyaw + unpitch
    for k=1:data.Nimages
        AHat = data.AHat(k,:) ;
        rawPsi(k)  = atan2(AHat(2),AHat(1)); % body angle with respect to x axis 
        rawBeta(k) = asin(AHat(3));  % body angle with respect to the horizon 
        rotM_YP(:,:,k) = eulerRotationMatrix(rawPsi(k), rawBeta(k),0) ; 
    end
   
end

% ---------------------------------
%% grab data from structure
AHat = data.AHat ;
bodyCM = data.bodyCM ;
spanHatR = data.rightSpanHats ;
spanHatL = data.leftSpanHats ;

Nimages = data.Nimages ;

% interpolate wing CM and tip?
if interpFlag
    [wingCMR, ~, wingTipR] = interpolateWingCM(data,'right') ;
    [wingCML, ~, wingTipL] = interpolateWingCM(data,'left') ;
else
    wingTipR = data.rightWingTips ;
    wingTipL = data.leftWingTips ;
    wingCMR = data.rightWingCM ;
    wingCML = data.leftWingCM ;
end

% -------------------------------------
%% transform data into body frame
% subtract off body CM
wingTipR = wingTipR - bodyCM ;
wingTipL = wingTipL - bodyCM ;

wingCMR = wingCMR - bodyCM ;
wingCML = wingCML - bodyCM ;

% rotate to align body axis with x hat
wingTipR_bodyFrame = nan(size(wingTipR)) ; 
wingTipL_bodyFrame = nan(size(wingTipL)) ; 
wingCMR_bodyFrame = nan(size(wingCMR)) ; 
wingCML_bodyFrame = nan(size(wingCML)) ; 

for i = 1:Nimages
    wingTipR_bodyFrame(i,:) = (squeeze(rotM_YP(:,:,i))*wingTipR(i,:)')' ; 
    wingTipL_bodyFrame(i,:) = (squeeze(rotM_YP(:,:,i))*wingTipL(i,:)')' ; 
    wingCMR_bodyFrame(i,:) = (squeeze(rotM_YP(:,:,i))*wingCMR(i,:)')' ; 
    wingCML_bodyFrame(i,:) = (squeeze(rotM_YP(:,:,i))*wingCML(i,:)')' ; 
end

% -----------------------------------------------------
%% get wing tip velocities in 
Ut_R = diff(wingTipR_bodyFrame) ;
Ut_R_xy = [Ut_R(:,1:2) zeros(length(Ut_R),1)] ;
Ut_R_xy = Ut_R_xy ./ repmat(myNorm(Ut_R_xy),1,3) ;
Ut_R_xy = [0, 0, 0; Ut_R_xy] ;

Ut_L = diff(wingTipL_bodyFrame) ;
Ut_L_xy = [Ut_L(:,1:2) zeros(length(Ut_L),1)] ;
Ut_L_xy = Ut_L_xy ./ repmat(myNorm(Ut_L_xy),1,3) ;
Ut_L_xy = [0, 0, 0; Ut_L_xy] ;

AHat_xy = [ones(Nimages,1), zeros(Nimages,2)] ; 

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
    rollVec = wingCML(goodRollPoints,:) - wingCMR(goodRollPoints,:) ;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
%% check for singularity in AHat to see if movie is likely a large pert
function singularityFlag = checkAHat(data)
    tol = 0.005 ;
    [minDist, ~] = min(abs(abs(data.AHat(:,3)) - 1)) ;
    singularityFlag = (minDist < tol) ;
end