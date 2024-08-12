% -------------------------------------------------------------------------
% function to generate a multipanel figure/image that will allow quick
% identification of movie type (i.e. pert type/direction/strength, when fly
% enters/leaves frame)
%
%{
rootPath = 'D:\Fly Data\Test\12_15102020\' ;
pathStruct = generatePathStruct(rootPath) ; 
MovNum = 10 ; 

[h_sum, ax_array] = makeMovieSummaryPanel(pathStruct, MovNum) ;
%}
% -------------------------------------------------------------------------
function [h_sum, ax_array] = makeMovieSummaryPanel(pathStruct, MovNum,...
    cineRangeFrames, pulseTiming)
% -------------------
%% inputs and params
if ~exist('cineRangeFrames', 'var') || isempty(cineRangeFrames)
    % if we don't know full cine range, try to get it from exprInfoStruct.
    exprInfoStruct = getExprInfo(pathStruct) ; 
    if isfield(exprInfoStruct, 'cineFrameStart') 
        cineRangeFrames = [exprInfoStruct.cineFrameStart, ...
            exprInfoStruct.cineFrameEnd] ; 
    else
        % faling that, just give a default value
        cineRangeFrames = [-723, 1086] ;  % full range of cine file (and thus mp4)
    end
end
if ~exist('pulseTiming', 'var') || isempty(pulseTiming)
    % similar to above, if we don't know pulseTiming, try to get it from
    % exprInfoStruct
    exprInfoStruct = getExprInfo(pathStruct) ; 
    if isfield(exprInfoStruct, 'optoPulseStart') 
       pulseTimingOpto = [ exprInfoStruct.optoPulseStart, ...
           exprInfoStruct.optoPulseEnd] ;
    else
        pulseTimingOpto = (1e-3).*[0, 50] ; 
    end
    if isfield(exprInfoStruct, 'magPulseStart') 
       pulseTimingMag = [ exprInfoStruct.magPulseStart, ...
           exprInfoStruct.magPulseEnd] ;
    else
        pulseTimingMag = (1e-3).*[15, 22] ; 
    end
    pulseTiming = [pulseTimingOpto ; pulseTimingMag] ; % in seconds -- default is for magnetic pulse
end

% constants for grabbing kinematic data
defineConstantsScript

% cine file extension (will be .cine for new cameras)
cine_file_ext = '.cine' ;

% order of camera views (from left to right) in mp4s
triptych_order = [2, 3, 1] ; % i.e. [xz, xy, yz]

% time interval between image samples
dt_samp = 0.015 ; % seconds

% read out experiment number from pathStruct
ExprNum = pathStruct.ExprNum ;

% params for pulse plots
pulseColors = lines(size(pulseTiming,1)) ; 
pulseAlpha = 0.2 ; 

% general figure size
figPosition = [460, 182, 1403, 800] ;
% --------------------------------------
%% load analysis data
% get path to data
suffixStr_in  = '_cleaned' ; % assume we're using cleaned data, but should be fine if not
[data, analysisOutput, ~, ~, errorFlag] = hierarchicalLoadData(pathStruct,...
    MovNum, [], suffixStr_in, true) ;

% if analysis file doesn't exist, exit
if errorFlag
    fprintf('Could not find path for movie %d \n', MovNum)
    h_sum = [] ;
    ax_array = [] ;
    return
end

% get start and stop time for movie; use to find t
startFrame = data.params.startTrackingTime ;  % first visible frame
endFrame = data.params.endTrackingTime ;
dt = (data.params.fps)^(-1) ;
t_start = dt*startFrame ;
t_end = dt*endFrame ;
t = t_start : dt : t_end ;

% use time for when we can track fly to get sampling times
t_samp1 = max([t_start - rem(t_start, dt_samp), -0.03]) ;
t_samp2 = min([t_end - rem(t_end, dt_samp), 0.05]) ;
t_samp = t_samp1 : dt_samp : t_samp2 ;

% read out body angles
anglesLabFrame = data.anglesLabFrame ;
angle_labels = {'body pitch', 'body roll'} ;
angle_inds = [THETAB, RHO] ;

% -------------------------------------------------------------------------
%% get mp4 movie file
% NB: this would probably be better with cine files so i'm going to add the
% skeleton of an option for using these, but will fill in later (don't
% often have cine files stored locally, so mp4 will likely be used more)
% -------------------------------------------------------------------------
% get paths to putative cine files
cine_path = pathStruct.root ;
cam_names = cellstr([data.params.cameraNames]);
cine_paths_full = cellfun(@(y) fullfile(cine_path, ...
    sprintf('%s_%03d%s',y, MovNum, cine_file_ext)), cam_names, ...
    'UniformOutput', false) ;

% check that cine files exist
cine_check = cellfun(@(y) exist(y, 'file'), cine_paths_full) ;
if all(cine_check)
    % --------------------------------------------------------------------
    % in this case, we have all the cine files. should be able to take
    % advantage of extant code to make this easier.
    % --------------------------------------------------------------------
    % clunky, but want to initialize storage for image structures
    im_struct_cell = cell(length(cam_names),1) ;
    
    % get root path for data type (input for "get_fly_snapshots")
    cine_path_split = strsplit(cine_path,filesep) ; 
    cine_root = strjoin(cine_path_split(1:end-1),filesep) ;
    
    % loop over cameras and use "get_fly_snapshots"
    for k = 1:length(cam_names)
        cam = cam_names{triptych_order(k)} ;
        im_struct = get_fly_snapshots(cine_root, ExprNum, MovNum, cam, ...
            t_samp, [], false, false, 1.0) ;
        
        % store image struct
        im_struct_cell{k} = im_struct ;
    end
    
else
    % --------------------------------------------------------------------
    % if we don't have cine files will need to use mp4 files (more likely
    % to be on local machine)
    % --------------------------------------------------------------------
    mp4_fn = sprintf('Expr_%d_movie_%03d.mp4', ExprNum, MovNum) ;
    mp4_path = pathStruct.mp4 ;
    mp4_path_full = fullfile(mp4_path, mp4_fn) ;
    
    % make sure mp4 file exists -- if not, that means no cine or mp4, so
    % quit
    if ~exist(mp4_path_full, 'file')
        fprintf('Error: could not find mp4 file %s \n', mp4_fn)
        h_sum = [] ;
        ax_array = [] ; 
        return
    end
    
    % if it does exist, loop over cameras and mp4 stills
    im_struct_cell = cell(length(cam_names),1) ;
    
    for k = 1:length(cam_names)
        cam = cam_names{triptych_order(k)} ;
        im_struct = get_fly_stills_from_mp4(pathStruct, MovNum, cam, ...
            t_samp, cineRangeFrames, analysisOutput) ;
        
        % store image struct
        im_struct_cell{k} = im_struct ;
    end
end

% ------------------------------------------
%% make sure we got images from all cameras
empty_check = cellfun(@(y) ~isfield(y,'image'), im_struct_cell) ; 
if any(empty_check)
    h_sum = [] ; 
    ax_array = [] ; 
    return
end
% --------------------------------------------------------------------
%% combine images obtained for each camera
% imSize = data.params.detectorLengthPix ;
im_comb_cell = cell(length(cam_names),1) ; 

% loop over cameras
for m = 1:length(cam_names)
    % current image struct
    im_struct = im_struct_cell{m} ;
    
    % current image size
    imSize = size(im_struct(1).image) ;
    
    % axlim for current camera
    axlim_curr = vertcat(im_struct.axlim) ; 
    
    % initialize sum image
    im_comb = uint8(zeros(imSize(1), imSize(2))) ;
    
    % loop over sampling times
    for n = 1:length(t_samp)
        % read out current image
        im_curr = im_struct(n).image ;
       
        % adjust contrast
        min_val = double(min(im_curr(:)))/255.0 ;
        im_adjust = imadjust(im_curr, [min_val, 1.0],[]) ;
        im_complement = imcomplement(im_adjust) ;
        
        % add processed image to summary image
        im_comb = imadd(im_comb, im_complement) ;
    end
    
    % at the end, take complement of image sum ...
    im_comb = imcomplement(im_comb) ;
    
    % ... and crop to tighter boundary
    indr1 = max([min(axlim_curr(:,3)), 1]) ;
    indr2 = min([max(axlim_curr(:,4)), imSize(1)]) ; 
    indr = indr1:indr2 ; 
    indc1 = max([min(axlim_curr(:,1)), 1]) ;
    indc2 = min([max(axlim_curr(:,2)), imSize(2)]) ; 
    indc = indc1:indc2 ;  
    im_comb = im_comb(indr, indc) ; 
    
    % store combined image in cell array
    im_comb_cell{m} = im_comb ; 
end

% -------------------------------------------------------------------------
%% put together summary panel
% initialize figure
h_sum = figure('PaperPositionMode', 'auto','OuterPosition', figPosition) ;

% subplot layout
n_rows = 3 ;
n_cols = 6 ;

% subplot gaps and margins
gap = [0.05, 0.05] ;
marg_h = [0.1, 0.1] ;
marg_w = [0.05, 0.05] ;

% storage for axes and plot counter
N_plots = 5 ;
ax_array = gobjects(N_plots,1) ;
cc = 1 ;

% ------------------------------------------
% first add image panels -- loop over cams
for k = 1:length(cam_names)
    % initialize subplot
    plt_rows = 1:2 ;
    plt_cols = (2*(k-1) + 1) : 2*k ;
    [r, c] = meshgrid(plt_rows, plt_cols) ;
    subplot_ind = sub2ind([n_cols, n_rows], c(:), r(:)) ;
    
    ax_array(cc) = subtightplot(n_rows, n_cols, subplot_ind, gap, ...
        marg_h, marg_w) ;
    set(ax_array(cc),'Parent', h_sum)
    
    % show image on axis
    image(ax_array(cc), im_comb_cell{k})  ;
    colormap(gray)
    
    % add title to image
    title(ax_array(cc), cam_names{triptych_order(k)})
    axis tight
    axis equal
    set(ax_array(cc), 'TickLength', [0,0], 'XTickLabel',[],'YTickLabel',[])
    
    % increment plot counter
    cc = cc + 1 ;
end

% -----------------------------------------------------
% then add body angle plots -- loop over angle indices
for m = 1:length(angle_inds)
    % read out curent angle and assoc label
    angleCurr = anglesLabFrame(:, angle_inds(m)) ;
    labelCurr = angle_labels{m} ;
    
    ylim = [min(angleCurr) - 5, max(angleCurr) + 5] ; 
    
    % initialize axis
    plt_rows = 3 ;
    plt_cols = (3*(m-1) + 1) : 3*m ;
    [r, c] = meshgrid(plt_rows, plt_cols) ;
    subplot_ind = sub2ind([n_cols, n_rows], c(:), r(:)) ;
    
    ax_array(cc) = subtightplot(n_rows, n_cols, subplot_ind, gap, ...
        marg_h, marg_w) ;
    set(ax_array(cc),'Parent', h_sum)
    hold(ax_array(cc), 'on') 
    
    % plot angles
    plot(ax_array(cc), t, angleCurr, 'k-') 
    
    % set/get axis limits
    set(ax_array(cc), 'xlim', [t(1), t(end)], 'ylim', ylim) 
    
    % --------------------------------------------
    % add markers for pulse timing
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    
    for q = 1:size(pulseTiming,1)
        % x axis corners of current pulse patch
        tsfvec = [pulseTiming(q,1), pulseTiming(q,2), pulseTiming(q,2), ...
            pulseTiming(q,1), pulseTiming(q,1)] ;
        
        % draw patch and set color properties
        hf = fill(ax_array(cc), tsfvec , avec,'y') ;
        set(hf,'facecolor',pulseColors(q,:),'facealpha',pulseAlpha,...
            'edgecolor','none','HandleVisibility','off') ;
    end
    
    % label axes
    xlabel(ax_array(cc), 'time (s)') 
    ylabel(ax_array(cc), [labelCurr ' (deg)'])
    title(ax_array(cc), labelCurr)
    
    % increment plot counter
    cc = cc + 1 ; 
end

% -------------------------------
%% add title to figure window
% create axis to house text
tit_ax = axes(h_sum, 'Position', [0.4, 0.93, 0.3, 0.06], 'Color', 'none') ;

% text for title
tit_str = strrep(sprintf('Expr_%d_mov_%03d', ExprNum, MovNum), '_',' ') ;
tit_txt = text(tit_ax, 0, 0, tit_str, 'Units','normalized', ...
    'VerticalAlignment','bottom','FontName', 'arial',...
    'FontSize', 24) ; 

% turn off axis rulers
axis(tit_ax, 'off')
end