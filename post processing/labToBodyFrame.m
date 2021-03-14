% -------------------------------------------------------------------------
% function to convert voxels, center of mass positions, and orientation
% vectors from the lab frame into the body frame of reference. 
%
% the body frame is defines as the long body axis pointing along x axis 
% with zero roll, pitched up by 45 degrees
% -------------------------------------------------------------------------
function data_out = labToBodyFrame(data_in, largePertFlag,thetaB0_rad, ...
    rollFlag) 
% ------------------------------------
%% params and inputs 
if ~exist('largePertFlag','var') || isempty(largePertFlag)
   largePertFlag = false ;  
end
if ~exist('thetaB0_rad','var') || isempty(thetaB0_rad)
    thetaB0_rad = -pi/4 ; % amount to pitch body up from x axis by (in strict body frame)
end
if ~exist('rollFlag','var') || isempty(rollFlag)
   rollFlag = true ; % if true, undo all rotations; if false, only unyaw and unpitch  
end
% initialize data_out structure
data_out = data_in ; 

% names of variables that require full transformation vs just rotation
trans_and_rot_names = {'rightWingCM', 'leftWingCM', 'rightWingTips',...
    'leftWingTips'} ; 
just_rot_names = {'rightChordHats', 'leftChordHats', 'rightSpanHats', ...
    'leftSpanHats' , 'AHat', 'rollVectors'} ; 

% get info for looping through data
N_frames = data_in.Nimages ; 
df = diff(data_in.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data_in.res,1)] ;
    
% -------------------------------------------------------------------
%% get transformations that will rotate/translate fly to body frame
% if these are already part of data, use those. otherwise, calculate
if isfield(data_out, 'rotM_YP') && isfield(data_out, 'rotM_roll')
    rotM_YP = data_out.rotM_YP ;
    rotM_roll = data_out.rotM_roll ;
else
    % estimate roll vectors if this hasn't been done already
    if ~isfield(data_out,'rollVectors')
        [rhoTimes, rollVectors] = estimateRollVector(data_out) ;
        data_out.rhoTimes = rhoTimes ;
        data_out.rollVectors = rollVectors ;
    end
    % get the set of matrices to unpitch + unyaw the fly, as well as the roll
    % angle
    [~, ~, ~, ~, ~, ~, ~, ~, ~, rotM_YP, rotM_roll, largePertFlag] =  ...
        calcAnglesRaw_Sam(data_out, false , largePertFlag) ;
end

% also define the rotation matrix to pitch the fly up by 45 degrees
rotM_theta0 = eulerRotationMatrix(0, thetaB0_rad, 0) ;

% loop through and create combined rotation matrix to unyaw, unpitch,
% unroll, and pitch up by 45 degrees
rotM_array = zeros(3,3,N_frames) ;
if rollFlag
    for k = 1:N_frames
        rotM_array(:,:,k) = rotM_theta0 * squeeze(rotM_roll(:,:,k)) * ...
            squeeze(rotM_YP(:,:,k)) ;
    end
else
    for k = 1:N_frames
        % version with no roll:
        rotM_array(:,:,k) = rotM_theta0 * squeeze(rotM_YP(:,:,k)) ;
    end
end

% -------------------------------------------------------------------
%% apply transformations to data to each frame
% NB: could probably use arrayfun here but w/e
for i = 1:N_frames
    % get rotation and translation for current frame
    rotM = squeeze(rotM_array(:,:,i)) ; 
    T = -1*data_in.bodyCM(i,:) ; 
    
    % first do the easy ones--just rotation
    for j = 1:length(just_rot_names)
       var_name = just_rot_names{j} ; 
       data_out.(var_name)(i,:) = (rotM*(data_out.(var_name)(i,:))')' ; % fucking transpose
    end
    
    % next do ones that require rotation AND translation
    for k = 1:length(trans_and_rot_names)
       var_name = trans_and_rot_names{k} ; 
       data_out.(var_name)(i,:) = (rotM*(data_out.(var_name)(i,:) + T)')' ; % fucking transpose
    end
    
    % finally, do the voxels, which aren't just one vector per frame
    row1 = frameStartInd(i) ; 
    row2 = frameEndInd(i) ; 
    data_out.res(row1:row2,2:4) = ...
        int16((rotM*double(data_out.res(row1:row2,2:4) + int16(T))')') ; % fucking double vs int16
end

% ------------------------------------------------------------
%% store the info that we'll need to get back to the lab frame
data_out.bodyCM_orig = data_in.bodyCM ; 
data_out.bodyFrameRotMats = rotM_array ; 
data_out.bodyCM = zeros(N_frames, 3) ; 
data_out.largePertFlag = largePertFlag ; 
%data_out.rollFlag = rollFlag ; 
%data_out.thetaB0_rad = thetaB0_rad ; 
end