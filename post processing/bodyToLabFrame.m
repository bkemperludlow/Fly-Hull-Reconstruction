% -------------------------------------------------------------------------
% function to convert voxels, center of mass positions, and orientation
% vectors from the body reference frame into the lab frame  
%
% the body frame is defines as the long body axis pointing along x axis 
% with zero roll, pitched up by 45 degrees
% -------------------------------------------------------------------------
function data_out = bodyToLabFrame(data_in)
% ----------------------------------
%% params
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

% initialize data_out struct
data_out = data_in ; 

% -----------------------------------------------------------------
%% make sure we have the necessary info for frame transformation
if ~isfield(data_in,'bodyFrameRotMats') || ~isfield(data_in,'bodyCM_orig')
   disp('Insufficient information to transform back to lab frame')
   return
end

bodyFrameRotMats = data_in.bodyFrameRotMats ; 
translationVecs = data_in.bodyCM_orig ; 

% ------------------------------------------------------------------
%% perform body->lab frame transformations
for i = 1:N_frames
    % get rotation matrix and translation vectors for current frame
    rotM = squeeze(bodyFrameRotMats(:,:,i))' ; % not the transpose/inverse
    T = translationVecs(i,:) ; 
   
    % first do the easy ones--just rotation
    for j = 1:length(just_rot_names)
       var_name = just_rot_names{j} ; 
       data_out.(var_name)(i,:) = (rotM*(data_out.(var_name)(i,:))')' ; % fucking transpose
    end
    
    % next do ones that require rotation AND translation
    for k = 1:length(trans_and_rot_names)
       var_name = trans_and_rot_names{k} ; 
       data_out.(var_name)(i,:) = (rotM*(data_out.(var_name)(i,:))')' + T ; % fucking transpose
    end
    
    % finally, do the voxels, which aren't just one vector per frame
    row1 = frameStartInd(i) ; 
    row2 = frameEndInd(i) ; 
    data_out.res(row1:row2,2:4) = ...
        int16((rotM*double(data_out.res(row1:row2,2:4))')' + T) ; % fucking double vs int16
    
end

% --------------------------------------------------------------
%% assign additional variables
data_out.bodyCM = translationVecs ; 
end