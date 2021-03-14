% -------------------------------------------------------------------------
% function to trim the time captured by a data structure
%
% time values in seconds!
%
% the motivation here is that some early fuck-ups with the angle can make
% the large perturbation calculations all messed up (due to wind-up error)
% -------------------------------------------------------------------------
function data_out = trimDataStruct(data_in, t_start_new, t_end_new, ...
    largePertFlag)
% -------------------------------
%% inputs 
if ~exist('t_start_new','var') || isempty(t_start_new)
    t_start_new = (1/data_in.params.fps)*data_in.params.startTrackingTime ; 
end
if ~exist('t_end_new','var') || isempty(t_end_new)
    t_end_new = (1/data_in.params.fps)*data_in.params.endTrackingTime ; 
end
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false ; 
end

% initialize data_out as a copy of data_in 
data_out = data_in ; 

% --------------------------------------------------------------------
%% lists of data structure fields that will need to be trimmed

% fields that contain data in the form 1xN (or Nx1)
scalarDataFields = {'rightChordTopProjections', 'leftChordTopProjections',...
    'diag11Right', 'diag12Right', 'diag21Left' , 'diag22Left' } ;  

% fields with data of the form NxM (where M > 1)
vectorDataFields = {'bodyCM', 'rightWingCM', 'leftWingCM' ,...
    'rightChordHats', 'leftChordHats', 'rightSpanHats',...
    'leftSpanHats', 'AHat', 'chord1AltHats', 'chord2AltHats', ...
    'rollVectors', 'bodyCM_old', 'AHat_old', 'rightWingTips',...
    'leftWingTips' };

% also NB: rather than deal with all the angle variables like phiL_amp,
% we'll just re-calculate the angles at the end

% other things we need to make sure to change:
%   -res and RESIDX (should be easy)
%   -data.params.startTrackingTime and data.params.endTrackingTime
%   -ignoreFrames
%   -Nimages
%   -rhoTimes
               
% --------------------------------------------------------------------
%% get time info from original data struct
frames_in = (data_in.params.startTrackingTime : data_in.params.endTrackingTime) ; 
t_in = (1/data_in.params.fps)*frames_in ; 

[~, t1_ind] = min(abs(t_in - t_start_new)) ; 
[~, t2_ind] = min(abs(t_in - t_end_new)) ; 

N_frames_new = length(t1_ind:t2_ind) ; 
N_frames_old = data_in.Nimages ; 
t1_ind_diff = t1_ind - 1 ; 
t2_ind_diff = N_frames_old - t2_ind + 1 ; 

frame1 = frames_in(t1_ind) ;
frame2 = frames_in(t2_ind) ; 
% -------------------------------------------------
%% trim data.params, data.Nimages, and ignoreFrames
fprintf('Trimming data... \n')

data_out.params.startTrackingTime = frame1 ; 
data_out.params.endTrackingTime = frame2  ; 
data_out.Nimages = N_frames_new ; 

if isfield(data_in,'ignoreFrames')
    ignoreFrames = data_in.ignoreFrames ; 
    if ~isempty(ignoreFrames)
       trim_idx = (ignoreFrames >= t1_ind) & (ignoreFrames <= t2_ind) ; 
       ignoreFrames = ignoreFrames(trim_idx) - t1_ind_diff ; 
    end
    data_out.ignoreFrames = ignoreFrames ; 
end

% -----------------------------------------------------
%% trim res and RESIDX
res = data_in.res ; 
RESIDX = data_in.RESIDX ;
res_trim_idx = (res(:,1) >= frame1) & (res(:,1) <= frame2) ; 

res = res(res_trim_idx,:) ; 
RESIDX = RESIDX(res_trim_idx,:) ; 

data_out.res = res ; 
data_out.RESIDX = RESIDX ; 

% ------------------------------------------------------
%% trim scalar and vector data structures
for i = 1:length(scalarDataFields)
   data_out.(scalarDataFields{i}) = ...
       data_in.(scalarDataFields{i})(t1_ind:t2_ind) ; 
end

for j = 1:length(vectorDataFields)
   data_out.(vectorDataFields{j}) = ...
       data_in.(vectorDataFields{j})(t1_ind:t2_ind,:) ; 
end

% -------------------------------------------------------
%% re-calculate angles to take care of the rest
% this is fucking stupid, but it's easier to save the data_out struct and
% then calculate angles using the main script (which takes a path as input)
% then do a work-around

fprintf('Re-calculating angles... \n')

% first estimate roll times and roll vectors (if not already there)
if isfield(data_in,'rhoTimes')
    rhoTimes = data_in.rhoTimes ;
    rho_trim_idx = (rhoTimes >= t1_ind) & (rhoTimes <= t2_ind) ;
    data_out.rhoTimes = rhoTimes(rho_trim_idx) - t1_ind_diff ;
    data_out.rollVectors = data_in.rollVectors(t1_ind:t2_ind,:) ;
else
    [rhoTimes, rollVectors] = estimateRollVector(data_out,largePertFlag, true) ;
    data_out.rhoTimes = rhoTimes ;
    data_out.rollVectors = rollVectors ;
end
            
% create temporary file
data_path = pwd ; 
data_filename = 'tttmp.mat' ;
temp_fn = fullfile(data_path, data_filename) ; 
save(temp_fn, 'data_out') ; 

% calculate angles
data_out = calcAnglesMain(temp_fn, largePertFlag, false, false) ; 

% delete temporary file
delete(temp_fn) ; 
end