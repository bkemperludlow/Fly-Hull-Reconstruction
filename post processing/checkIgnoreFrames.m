% -------------------------------------------------------------------------
% function to check if, after doing some post processing, we need to remove
% some frames from the "ignoreFrames" list -- e.g., if we correct roll and
% do some L/R swapping of wings, maybe some of the frames fit well now.
%
% NB: this would be easier if i had stored the reason for ignoring a given
% frame (e.g. bad clustering vs an attempted L/R swap that didn't work out.
% should maybe revise that in the future.
%
% going super basic for now (might have a better idea later, but this is
% all i have for not. basically just look at wing cm  -- if
% frames now seem okay on that front, remove them from ignore frames.
% ideally this should be pretty conservative -- we don't want to
% unnecessarily remove ignoreFrames
% -------------------------------------------------------------------------
function data = checkIgnoreFrames(data)
% ---------------------------------------
% check wing centers of mass
[~, good_idx_R, ~] = interpolateWingCM(data,'right') ;
[~, good_idx_L, ~] = interpolateWingCM(data,'left') ;

good_cm_ind = find(good_idx_R & good_idx_L) ;

% read out current ignore frames
ignoreFrames = data.ignoreFrames ; 

% remove ignore frames that are also found in good_cm_ind
ignoreFramesNew = setdiff(ignoreFrames, good_cm_ind) ; 

% reassign value of ignoreFrames 
data.ignoreFrames = ignoreFramesNew ; 

end