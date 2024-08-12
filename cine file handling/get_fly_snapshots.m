%--------------------------------------------------------------------------
% function to grab images from a cine file at specified times. Called in
% "plot_fly_snapshots.m" to create image that combines multiple stills from
% a single video
%
% Example usage:
%   im_struct = get_fly_snapshots('D:\Fly Data\Opto Silencing\', 4, 139, 'XY', ...
%       [-0.01, 0, 0.01, 0.02], [], false, true, 1.0)

%{
im_struct = get_fly_snapshots('D:\Box Sync Old\Haltere Experiment\', 22, 18, 'XY', ...
    [-0.01, 0, 0.01, 0.02], [], false, true, 1.0)
%}

%--------------------------------------------------------------------------
function im_struct = ...
    get_fly_snapshots(cinRoot, ExprNum, MovNum, cam, snapShotTimes, ...
    savePath, saveFlag, trimFlag, scale, fileExt)
% -------------------------------------------------------------------------
%% params and inputs
if ~exist('savePath','var') || isempty(savePath)
    savePath = [] ;
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = false ;
end
if ~exist('trimFlag', 'var') || isempty(trimFlag)
    trimFlag = true ;
end
if ~exist('scale','var') || isempty(scale)
    scale = 1.0 ; 
end
if ~exist('fileExt','var') || isempty(fileExt)
    fileExt = '.cine' ; % '.cin' | '.cine'
end

im_struct = struct() ;
DELTA = 18 ;

% -------------------------------------------------------------------------
%% get indices for cameras
if strcmpi(cam, 'XY')
    camNum = 3 ;
elseif strcmpi(cam, 'XZ')
    camNum = 2 ;
elseif strcmpi(cam, 'YZ')
    camNum = 1 ;
else
    keyboard
end

cam_name_cell = {'yz','xz','xy'} ;

% -------------------------------------------------------------------------
%% locate cine file and associated analysis
% first find folder corresponding to experiment number
exprDir = dir(cinRoot) ; 
exprDir = exprDir(3:end) ;

exprFolderInd = arrayfun(@(x) strcmp(x.name(1:2), ...
    num2str(ExprNum,'%02d')), exprDir) ;
exprFolder = fullfile(exprDir(exprFolderInd).folder,...
    exprDir(exprFolderInd).name) ;

% then get cine file path
cinFilename = fullfile(exprFolder, [cam_name_cell{camNum} '_' ...
    num2str(MovNum,'%03d') fileExt]) ;
%metaData = getCinMetaData(cinFilename) ;
cindata  = myOpenCinFile(cinFilename) ;

% then generate a pathStruct to find analysis info
pathStruct = generatePathStruct(exprFolder) ; 
analysisPath = findMovAnalysisPath(pathStruct, MovNum) ;
analysisPathFull = fullfile(analysisPath, ['Expr_' num2str(ExprNum) '_mov_' ...
    num2str(MovNum,'%03d') '_results.mat']) ;

analysisPathRecon =  fullfile(analysisPath, 'hullRecon.mat') ;
% -------------------------------------------------------------------------
%% load analysis data and get relevant info
if exist(analysisPathFull,'file')
    analysis_in = load(analysisPathFull) ;
    
    t_start = analysis_in.params.startTrackingTime ;
    t_end = analysis_in.params.endTrackingTime ;
    frames = (t_start - DELTA) : (t_end + DELTA) ;
    time_sec = (1/8000)*frames ;
    
    allAxlim = analysis_in.(['allAxlim_' cam_name_cell{camNum}]) ;
    if isfield(analysis_in, 'allBG')
        bg = squeeze(analysis_in.allBG(camNum,:,:)) ;
    elseif isfield(analysis_in, 'allBGcell')
        bg = analysis_in.allBGcell{camNum} ; 
    else
        fprintf('Error: cannot find background in %s \n', analysisPathFull)
        keyboard
    end
    
elseif ~exist(analysisPathFull,'file') && exist(analysisPathRecon,'file')
    hullRecon = load(analysisPathRecon) ; 
    binaryThresh = load(fullfile(analysisPath,'binaryThresh.mat')) ; 
    BG = load(fullfile(analysisPath,'BG.mat')) ; 
    
    t_start = hullRecon.params.startTrackingTime ;
    t_end = hullRecon.params.endTrackingTime ;
    frames = (t_start - DELTA) : (t_end + DELTA) ;
    time_sec = (1/8000)*frames ;
    
    allAxlim = binaryThresh.(['allAxlim_' cam_name_cell{camNum}]) ;
    bg = squeeze(BG.allBG(camNum,:,:)) ;
else
    disp('error loading analysis results')
    return
end
% -------------------------------------------------------------------------
%% loop through frames and grab images, process, then add to im_struct
for k = 1:length(snapShotTimes)
    % get proper time
    t_sec_curr = snapShotTimes(k) ;
    [~,t_ind] = min(abs(time_sec - t_sec_curr)) ;
    
    % read in raw image and subtract background
    im1 = myReadCinImage(cindata, frames(t_ind)) ;
    im2 = imsubtract(bg, im1) ;
    
    %frame_ind = find(time_sec == t_sec_curr,1,'first') ;
    % trim to axis limits?
    if trimFlag
        axlim_curr = round(allAxlim(t_ind,:)) ;
        
        im3 = im2(axlim_curr(3):axlim_curr(4),...
            axlim_curr(1):axlim_curr(2)) ;
    else
        im3 = im2 ;
    end
    
    % use a multi-level threshold to remove any background noise
    levels = multithresh(im3,3) ; 
    temp_im = imquantize(im3,levels) ; 
    ind_low = (temp_im == 1) ; 
    
    im4 = im3 ; 
    im4(ind_low) = 0 ; 
    im5 = imcomplement(im4) ; 
    
    % scale image up to reduce pixelation?
    if scale ~= 1.0 
        im5 = imresize(im5, scale) ; 
        im5 = imgaussfilt(im5, 1.0) ; 
    end
    
    % add image to structure along with time info
    im_struct(k).image = im5 ;
    im_struct(k).tsec = t_sec_curr ;
    im_struct(k).tframe = frames(t_ind) ;
    im_struct(k).frame = t_ind ;
    im_struct(k).axlim = round(allAxlim(t_ind,:)) ;
    
    % save image?
    if saveFlag
        fn = fullfile(savePath, ['Expr_' num2str(ExprNum) '_mov_' ...
            num2str(MovNum,'%03d') '_' cam '_' num2str(t_ind,'%04d')...
            '.png']) ; 
        imwrite(im5, fn)
    end
    
end

%% close cine file
myCloseCinFile(cindata) ;
end