function [bg, tin, tout, xcm, ycm] = findBG_MOG_mk2(cinFilename,minPixArea)
% find the background image for the given movie using a mixture of
% gaussians model. Draws heavily for multi object tracking example
%
% minPixArea is the minimum blob size to qualify for
% detection. default value of 20 seems to work okay
%
% bg is the image, 
% tin and tout are estimations for the time the fly goes in and out the
% frame.
% xcm and ycm are a ROUGH estimate for the center-of-mass of the fly

phantomSDK_setPath ;
LoadPhantomLibraries;

warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString')

DEBUG_FLAG = false;

%------------------------------
% params

if (~exist('minPixArea','var'))
    minPixArea = 20 ; 
end
% tracker
assignmentThresh = 30 ; 
numCoastingUpdates = 10 ;
confirmationParams = [9, 10] ; 
%------------------------------
% get movie data
metaData = getCinMetaData(cinFilename) ;
cindata  = myOpenCinFile(cinFilename) ;

% movie characteristics 
imWidth = metaData.width ;
imHeight = metaData.width ;
imMidPoint_x = round(imWidth/2) ; 
imMidPoint_y = round(imHeight/2) ; 

N    = double(metaData.lastImage) - double(metaData.firstImage) + 1 ;
%df   = zeros(N,1) ;
minTrackLength = 200 ; 
tvec = (double(metaData.firstImage)) : double(metaData.lastImage) ;
%SE = strel('disk',4) ; % se size was 2

%----------------------------------
% initialize arrays 

tracks_struct = struct() ; 
t0_ind = find(tvec == 0) ; 

%----------------------------------

if (~DEBUG_FLAG)
    hbar = waitbar(0,['findBG: Processing images in ' cinFilename(end-9:end)]) ;
end

%--------------------------------------------------------------------------
% set up tracking objects
detectorObjects = setupDetectorObjects(minPixArea) ; 
tracker = multiObjectTracker(...
    'FilterInitializationFcn', @flyTrackFilter, ...
    'AssignmentThreshold', assignmentThresh, ...
    'NumCoastingUpdates', numCoastingUpdates, ...
    'ConfirmationParameters', confirmationParams ...
    );

if DEBUG_FLAG 
    h_debug = figure ; 
    ax = gca ; 
end

%--------------------------------------------------------------------------
% loop through video with mixture of gaussians detector

c = 1 ;

for it = metaData.firstImage : metaData.lastImage 
    frame = myReadCinImage(cindata, it) ;
    
    [detections, mask] = detectObjects(detectorObjects, frame, c, minPixArea);
    confirmedTracks = updateTracks(tracker, detections, c);
    tracks_struct(c).confirmedTracks = confirmedTracks ; 
    
    if DEBUG_FLAG
        displayTrackingResults(ax, confirmedTracks, frame, mask,it);
    end
    c = c+ 1 ; 
    if (~DEBUG_FLAG)
        waitbar(c/N, hbar) ;
    end
end

if (~DEBUG_FLAG)
    close(hbar)
end

%----------------------------------
% separate out different tracks
empty_tracks_ind = arrayfun(@(x) isempty(x.confirmedTracks), tracks_struct) ; 
tracks_struct_nonempty = tracks_struct(~empty_tracks_ind) ; 

trackIDs = arrayfun(@(x) [x.confirmedTracks.TrackID], tracks_struct_nonempty,...
    'UniformOutput',false) ;
trackIDs_unique = unique(cell2mat(trackIDs)) ; 

xcm_mat = nan(N, length(trackIDs_unique)) ; 
ycm_mat = nan(N, length(trackIDs_unique)) ; 
area_mat = nan(N, length(trackIDs_unique)) ; 

for j = 1:length(tracks_struct_nonempty)
    confirmedTracks = tracks_struct_nonempty(j).confirmedTracks ; 
    for k = 1:length([confirmedTracks.TrackID])
       trackID_curr = confirmedTracks(k).TrackID ;
       state_curr = confirmedTracks(k).State ;
       area_curr = confirmedTracks(k).ObjectAttributes{1}{2} ; 
       frame_num_curr = confirmedTracks(k).Time ; 
       
       mat_ind = find(trackIDs_unique == trackID_curr) ; 
       
       xcm_mat(frame_num_curr,mat_ind) = state_curr(1) ; 
       ycm_mat(frame_num_curr,mat_ind) = state_curr(3) ; 
       area_mat(frame_num_curr,mat_ind) = area_curr ; 
    end
end

% filter out trajectories that are too short
track_length = sum(~isnan(xcm_mat),1) ; 

xcm_mat = xcm_mat(:,(track_length > minTrackLength)) ; 
ycm_mat = ycm_mat(:,(track_length > minTrackLength)) ; 
area_mat = area_mat(:,(track_length > minTrackLength)) ; 
%----------------------------------
% try to find "main" fly track
diff_from_center = sqrt((xcm_mat(t0_ind,:) - imMidPoint_x).^2 + ...
    (ycm_mat(t0_ind,:) - imMidPoint_y).^2) ;  

[~, mainFlyInd] = min(diff_from_center) ; 

xcm = xcm_mat(:,mainFlyInd) ; 
ycm = ycm_mat(:,mainFlyInd) ; 
area_fly = area_mat(:,mainFlyInd) ; 
%---------------------------------
out = true(N,1) ; 
for q = 1:size(xcm_mat,2)
   out = out & isnan(xcm_mat(:,q)) & isnan(ycm_mat(:,q)) ;  
end
%out = (xcm == 0) & (ycm == 0) ; 
%entry_ind = find(diff(out) == -1) + 1 ; 
%entry_ind = entry_ind(entry_ind < round(N/2)) ; 

% decide when the fly was in FOV
%tin = tvec(find(~out,1,'first')) ; %this could be changed so that it goes from zero outwards
%tout = tvec(find(~out,1,'last')) ;
%-------------------------------------

% if there is a time in which the fly is outside the FOV, take that as
% background
%out = medfilt1(double(out),10) ;
bgIndex = find(out,1,'first') ; % find the first "out" image and make it the BG image
firstEntranceIdx = find(~out,1,'first') ; 

if (~isempty(bgIndex)) && (firstEntranceIdx > (assignmentThresh + 1))
    %bg = ReadCineFileImage(cinFilename, tvec(bgIndex), false);
    bg  = myReadCinImage(cindata, tvec(bgIndex)) ;
else
     % if the fly is always in FOV, fabricate a BG image (automatically, yay)
   
    % find the two most distant points along the trajectory
    D = squareform(pdist([xcm ycm])) ;
    [m1, ~]   = max(D) ;
    [~, ind2] = max(m1) ;
    [~, ind3] = max(D(:, ind2)) ;
    % max is in D(ind3, ind2) ;
    clear D
    
    t1 = ind3;
    t2 = ind2;
    
    xmid = round ( ( xcm(t1) + xcm(t2) ) /2 );
    ymid = round ( ( ycm(t1) + ycm(t2) ) /2 );
    
    % find if the larger difference is along x or y
    dx = abs(xcm(t1) - xcm(t2)) ;
    dy = abs(ycm(t1) - ycm(t2)) ;
    
    % read the two images
    %imt1 = ReadCineFileImage(cinFilename, tvec(t1), false);
    %imt2 = ReadCineFileImage(cinFilename, tvec(t2), false);
    
    imt1 = myReadCinImage(cindata, tvec(t1)) ;
    imt2 = myReadCinImage(cindata, tvec(t2)) ;
    
    bg = imt1 ;
    
    if (dx>=dy)
        if (xcm(t1)<xmid) % if fly is left to the mid point at t1, take the left part from t2
            bg(:, 1:xmid) = imt2(:,1:xmid) ;
        else % if fly is right to the midpoint at t1, take the right part from t2
            bg(:, xmid:end) = imt2(:, xmid:end) ;
        end
    else
        if (ycm(t1)<ymid) % if fly is above to the mid point at t1, take the top part from t2
            bg(1:ymid,:) = imt2(1:ymid,:) ;
        else % if the fly is below the midpoint at t1, take the bottom part from t2
            bg(ymid:end,:) = imt2(ymid:end,:) ;
        end
    end
end
%------------------------------------------------------------------------
% find first and last images that contain THE FULL FLY (or most of it, at
% least)
flyPixAreaGuess = median(area_fly(~isnan(area_fly))) ; 

full_fly_frames = false(N,1) ; 
full_fly_frames(area_fly >= flyPixAreaGuess) = true ; 

tin_ind = find(diff(full_fly_frames) == 1,1,'first') + 1 ;
tout_ind = find(diff(full_fly_frames) == -1,1,'last')  ; 

if isempty(tin_ind)
    tin_ind = 1 + 80 ;
end
if isempty(tout_ind)
    tout_ind = length(tvec) - 80 ;
end

tin = tvec(tin_ind) ;
tout = tvec(tout_ind) ;

% inFrame_ind = false(size(tvec)) ; 
% inFrame_ind(tin_ind:tout_ind) = true ; 
% zeroCM_ind = (xcm == 0) & (ycm == 0) ; 
% zeroCM_ind = zeroCM_ind' ; 
% 
% toInterp_ind = inFrame_ind & zeroCM_ind ; 
% toFit_ind = inFrame_ind & ~zeroCM_ind ; 
% 
% xcm(toInterp_ind) = interp1(tvec(toFit_ind), xcm(toFit_ind),tvec(toInterp_ind)) ; 
% ycm(toInterp_ind) = interp1(tvec(toFit_ind), ycm(toFit_ind),tvec(toInterp_ind)) ; 

%------------------------------------------------------------------------
if (DEBUG_FLAG)
    figure ;
    plot(tvec, df,'k.-') ;
    hold on ;
    plot(tvec, out*50,'r-','linewidth',2) ;
    %plot(tvec, in*50,'g-','linewidth',2) ;
    hold off ;
    
    figure ;
    plot(xcm, ycm,'bo-') ;
    axis equal ;
    
    %keyboard ;
end

myCloseCinFile(cindata) ;

warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString')

return 
end

%==========================================================================

function detectorObjects = setupDetectorObjects(minBlobArea)
% Create System objects for foreground detection and blob analysis

% The foreground detector segments moving objects from the
% background. It outputs a binary mask, where the pixel value of 1
% corresponds to the foreground and the value of 0 corresponds to
% the background.

detectorObjects.detector = vision.ForegroundDetector('NumGaussians', 5, ...
    'NumTrainingFrames', 50, 'MinimumBackgroundRatio', 0.7,...
    'InitialVariance',60^2);

% Connected groups of foreground pixels are likely to correspond to
% moving objects.  The blob analysis System object finds such
% groups (called 'blobs' or 'connected components') and computes
% their characteristics, such as their areas, centroids, and the
% bounding boxes.

detectorObjects.blobAnalyzer = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...
    'AreaOutputPort', true, 'CentroidOutputPort', true, ...
    'MinimumBlobArea', minBlobArea);
end

%==========================================================================

function [detections, mask_filt2] = detectObjects(detectorObjects, frame,...
    frameCount,minPixArea)
% Expected uncertainty (noise) for the blob centroid.
measurementNoise = 20*eye(2);
%minPixArea = 20 ; 
% Detect foreground.
mask = detectorObjects.detector.step(frame);

% Apply morphological operations to remove noise and fill in holes.
mask_filt = bwareaopen(mask, minPixArea) ; 
mask_filt2 = imopen(mask_filt,strel('disk',2)) ; 
mask_filt2 = imclose(mask_filt2,strel('disk',10)) ; 
mask_filt2 = imfill(mask_filt2,'holes') ; 
    
% Perform blob analysis to find connected components.
[areas, centroids, bboxes] = detectorObjects.blobAnalyzer.step(mask_filt2);

% Formulate the detections as a list of objectDetection objects.
numDetections = size(centroids, 1);
detections = cell(numDetections, 1);
for i = 1:numDetections
    detections{i} = objectDetection(frameCount, centroids(i,:), ...
        'MeasurementNoise', measurementNoise, ...
        'ObjectAttributes', {bboxes(i,:),areas(i)});
end
end

%==========================================================================

function filter = flyTrackFilter(detection)
% Initialize a Kalman filter for this example.

% Define the initial state.
state = [detection.Measurement(1); 0; detection.Measurement(2); 0];

% Define the initial state covariance.
stateCov = diag(0.5*[50, 50, 50, 50]);

% Create the tracking filter.
filter = trackingKF('MotionModel', '2D Constant Velocity', ...
    'State', state, ...
    'StateCovariance', stateCov, ...
    'MeasurementNoise', detection.MeasurementNoise(1:2,1:2) ...
    );
end

%==========================================================================

function displayTrackingResults(ax_handle, confirmedTracks, frame,...
    mask, frameTime)
% Convert the frame and the mask to uint8 RGB.
frame = im2uint8(frame);
mask = uint8(repmat(mask, [1, 1, 3])) .* 255;

if ~isempty(confirmedTracks)
    % Display the objects. If an object has not been detected
    % in this frame, display its predicted bounding box.
    numRelTr = numel(confirmedTracks);
    boxes = zeros(numRelTr, 4);
    ids = zeros(numRelTr, 1, 'int32');
    predictedTrackInds = zeros(numRelTr, 1);
    for tr = 1:numRelTr
        % Get bounding boxes.
        boxes(tr, :) = confirmedTracks(tr).ObjectAttributes{1}{1};
        
        % Get IDs.
        ids(tr) = confirmedTracks(tr).TrackID;
        
        if confirmedTracks(tr).IsCoasted
            predictedTrackInds(tr) = tr;
        end
    end
    
    predictedTrackInds = predictedTrackInds(predictedTrackInds > 0);
    
    % Create labels for objects that display the predicted rather
    % than the actual location.
    labels = cellstr(int2str(ids));
    
    isPredicted = cell(size(labels));
    isPredicted(predictedTrackInds) = {' predicted'};
    labels = strcat(labels, isPredicted);
    
    % Draw the objects on the frame.
    frame = insertObjectAnnotation(frame, 'rectangle', boxes, labels);
    
    % Draw the objects on the mask.
    mask = insertObjectAnnotation(mask, 'rectangle', boxes, labels);
end

% Display the mask and the frame.
imshowpair(frame,mask,'montage','Parent',ax_handle) ; 
title(num2str(frameTime))
pause(1/50)
end

%==========================================================================
