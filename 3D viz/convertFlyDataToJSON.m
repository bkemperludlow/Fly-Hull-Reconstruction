% -------------------------------------------------------------------------
% function to convert fly data structure to .JSON format for compatibility 
% with Maya and potentially other 3D animation software
%
% because we can't just convert the whole data structure, need to pick
% which things we'd like to convert. For animation, it seems like the most
% relevant things are:
%   1) body CM 
%   2) body orientation rotation matrix
%   3) L/R wing rotation matrices
%   4) Time
%   5) wing CM? probably not necessary, since, for the animation, we assume
%   a fixed distance
%
% -------------------------------------------------------------------------
function [bodyCM_json, bodyRot_json, wingRotR_json, wingRotL_json] = ...
    convertFlyDataToJSON(data, savePath, largePertFlag)
% ------------------------------
%% inputs
if ~exist('savePath','var') || isempty(savePath)
   savePath = 'D:\Fly Data\For Vision Research\' ; 
end
if ~exist('largePertFlag','var') || isempty(largePertFlag)
   largePertFlag = guessLargePert(data) ; 
end
debugFlag = false ; % make vector plot to check rotation matrices?
plotFlag = false ; % just for angle calculations
theta0 = pi/4 ; 

% -----------------------------
%% read relevant data from structure
% params
voxelSize = data.params.voxelSize ; 
Nframes = data.Nimages ; 

% time info
tFrames = data.params.startTrackingTime : data.params.endTrackingTime ; 
tsec = (1/8000)*tFrames ; 

% body center of mass
bodyCM_vox = data.bodyCM ; 
bodyCM = voxelSize*bodyCM_vox ; 

% wing angles (body frame)
[~, smooth_anglesMat_R, ~, ~, ~] = smoothWingAngles(data, 'R') ;
[~, smooth_anglesMat_L, ~, ~, ~] = smoothWingAngles(data, 'L') ; 
% convet to radians
anglesMat_R = (pi/180) * smooth_anglesMat_R ; 
anglesMat_L = (pi/180) * smooth_anglesMat_L ; 

% rotation matrix to get fly to pitch = 45 degrees
rotM_theta0 = eulerRotationMatrix(0, -theta0, 0) ; 

% ----------------------------------------------------------
%% get wing angle rotation matrices
rotM_wingR = zeros(3,3,Nframes) ; 
rotM_wingL = zeros(3,3,Nframes) ; 
for i = 1:Nframes
    % get rotation matrix based purely on euler angles
    tempR = eulerRotationMatrix(-1*anglesMat_R(1,i), anglesMat_R(2,i), ...
        anglesMat_R(3,i)) ; 
    tempL = eulerRotationMatrix(anglesMat_L(1,i), anglesMat_L(2,i), ...
        anglesMat_L(3,i)) ; 
    
    % include rotation to stroke plane
    rotM_wingR(:,:,i) =  rotM_theta0' * tempR'  ;
    rotM_wingL(:,:,i) =  rotM_theta0' * tempL'  ;
end
% ----------------------------------------------------------------
%% (re) calculate BODY angle rotation matrices 
[~, ~, ~, ~,~, ~, ~, ~, ~, rotM_YP, rotM_roll] = ...
    calcAnglesRaw_Sam(data, plotFlag ,largePertFlag) ; 
rotM_body = zeros(size(rotM_YP)) ; 
for m = 1:size(rotM_YP,3)
    rotM_body(:,:,m) = (rotM_roll(:,:,m) * rotM_YP(:,:,m))' ; 
end

% ----------------------------------------------------------
%% test rotation matrices?
if debugFlag
    ind = 50 ; 
   R_wingR = squeeze(rotM_wingR(:,:, ind)) ; 
   R_wingL = squeeze(rotM_wingL(:,:, ind)) ;
   R_body = squeeze(rotM_body(:,:,ind)) ; 
   
   x_hat = [1; 0; 0] ; 
   bodyAxisHat = R_body*x_hat ; 
   spanHatR = R_body*R_wingR*x_hat ; 
   spanHatL = R_body*R_wingL*x_hat ; 
   figure ;
   hold on
   plot3([0, bodyAxisHat(1)] + bodyCM(ind, 1), ...
       [0, bodyAxisHat(2)] + bodyCM(ind, 2), ...
        [0, bodyAxisHat(3)] + bodyCM(ind, 3), ...
        'k-', 'LineWidth', 2)
    plot3([0, spanHatR(1)] + bodyCM(ind, 1) + 0.5*bodyAxisHat(1), ...
       [0, spanHatR(2)] + bodyCM(ind, 2) + 0.5*bodyAxisHat(2), ...
        [0, spanHatR(3)] + bodyCM(ind, 3) + 0.5*bodyAxisHat(3), ...
        'r-', 'LineWidth', 2)
    plot3([0, spanHatL(1)] + bodyCM(ind, 1) + 0.5*bodyAxisHat(1), ...
       [0, spanHatL(2)] + bodyCM(ind, 2) + 0.5*bodyAxisHat(2), ...
        [0, spanHatL(3)] + bodyCM(ind, 3) + 0.5*bodyAxisHat(3), ...
        'b-', 'LineWidth', 2)
    plot3(bodyAxisHat(1) + bodyCM(ind, 1), ...
        bodyAxisHat(2) + bodyCM(ind, 2), ...
        bodyAxisHat(3) + bodyCM(ind, 3), 'ko','MarkerFaceColor','g')
    axis equal 
    grid on 
    box on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    keyboard
end

% --------------------------------------------------------------------
%% convert to JSON and save
bodyCM_withT = [tsec' , bodyCM] ; 
bodyCM_json = savejson('bodyCM', bodyCM_withT,'NestArray',0) ;
save(fullfile(savePath, 'bodyCM.json'),'bodyCM_json')

bodyRot_json = savejson('bodyRot', rotM_body,'NestArray',0) ;
save(fullfile(savePath, 'bodyRot.json'),'bodyRot_json')

wingRotR_json = savejson('wingRotR', rotM_wingR,'NestArray',0) ;
save(fullfile(savePath, 'wingRotR.json'),'wingRotR_json')

wingRotL_json = savejson('wingRotL', rotM_wingL,'NestArray',0) ;
save(fullfile(savePath, 'wingRotL.json'),'wingRotL_json')

end