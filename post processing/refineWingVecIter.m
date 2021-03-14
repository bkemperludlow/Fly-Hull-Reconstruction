% -------------------------------------------------------------------------
% function to correct errors in span (and chord?) estimation by looking at
% the time series of vector orientations. May need to use images/voxels as
% well
% -------------------------------------------------------------------------
function data_out = refineWingVecIter(data_in, interpType, bodyFrameFlag)
% ------------------------
%% params and inputs
if ~exist('interpType','var') || isempty(interpType) 
    interpType = 'spline' ; % 'spline' or 'rot'
end
if ~exist('bodyFrameFlag','var') || isempty( bodyFrameFlag) 
     bodyFrameFlag = false ; % check wing vectors in body frame of reference
end

% minimum value for dot product of spans in successive frames
dotThreshSpan = 0.90 ; 
dotThreshChord = 0.00 ; 
% method of interpolation
% max number of iterations to converge on good set of wing vecs
N_iter_max = 50 ; 
% number of frames
N_frames = data_in.Nimages ;
% indicators of wing side to loop through
wing_strs = {'right','left'} ;  %{'right', 'left'} ; 
wing_num_strs = {'1', '2'} ;  
% unit vectors
x_hat = [1, 0, 0] ; 
% ------------------------------------------------------
%% initialize output
if bodyFrameFlag
    largePertFlag = guessLargePert(data_in) ; 
    data_out = labToBodyFrame(data_in, largePertFlag) ; 
else
    data_out = data_in ; 
end
ignoreFrames_new = [] ; 

% also get any ignoreFrames from data struct, if they exist
if isfield(data_in, 'ignoreFrames')
    ignoreFrames = data_in.ignoreFrames ; 
else
    ignoreFrames = [] ; 
end
% ------------------------------------------------------
%% get wing vectors and compare adjacent frames
% loop over left and right side, since procedure should be the same
for ws = 1:length(wing_strs)
    % read in wing vectors from data structure
    spanHats = data_out.([wing_strs{ws} 'SpanHats']) ;
    chordHats = data_out.([wing_strs{ws} 'ChordHats']) ;
    chordAltHats = data_out.(['chord' wing_num_strs{ws} 'AltHats']) ;
    
    fprintf('Correcting spans for %s wing... \n', wing_strs{ws})
    % get rotation matrices for each frame
    rotM_span_array = zeros(3,3,N_frames) ;
    span_nan_idx = any(isnan(spanHats(2:end,:)),2) ; 
    loop_val = find(~span_nan_idx)' + 1 ; 
    for i = loop_val
        spanVR = vrrotvec(x_hat, spanHats(i,:)) ; 
        rotM_span_array(:,:,i) = vrrotvec2mat(spanVR) ; 
    end
    
    % ---------------------------------------------
    %% first deal with spans
    % find frames where dot product between adjacent frames is off and
    % interpolate through them. do this until all frames look good or we
    % reach iteration limit
    [spanHats, ignoreFramesSpan] = iterInterpWingVecs(spanHats, ...
        dotThreshSpan, rotM_span_array, N_iter_max, ignoreFrames,interpType) ;
    data_out.([wing_strs{ws} 'SpanHats']) = spanHats ; 
    
    fprintf('Span correction complete for %s wing \n', wing_strs{ws})
    % --------------------------------------------------------
    %% next check chords
    % first see if alternate chords give a better fit in any frames
%     chordDot = dot(chordHats(1:(end-1),:), chordHats(2:end,:), 2) ;
%     chordAltDot = dot(chordHats(1:(end-1),:), chordAltHats(2:end,:), 2) ;
%     
%     swap_idx = find(chordAltDot > chordDot) ;
%     %disp(sum(1-chordDot)) 
%     cc = 0 ; 
%     while (numel(swap_idx) > 0) && (cc < N_iter_max) 
%         cc = cc + 1 ; 
%         tmp = chordHats ; 
%         chordHats(swap_idx + 1, :) = chordAltHats(swap_idx + 1,:) ; 
%         chordAltHats(swap_idx + 1,:) = tmp(swap_idx + 1,:) ; 
%         
%         chordDot = dot(chordHats(1:(end-1),:), chordHats(2:end,:), 2) ;
%         chordAltDot = dot(chordHats(1:(end-1),:), chordAltHats(2:end,:), 2) ;
%         
%         swap_idx = find(chordAltDot > chordDot) ;
%         
%     end
    %disp(sum(1-chordDot)) 
    
    %-------------------------------------------------------------------
    %% then use same trick with span to interpolate through bad points
    % get (chord) matrices for each frame
    fprintf('Correcting chords for %s wing... \n', wing_strs{ws})
    rotM_chord_array = zeros(3,3,N_frames) ;
    chord_nan_idx = any(isnan(chordHats(2:end,:)),2) ; 
    loop_val = find(~chord_nan_idx)' + 1 ; 
    for i = loop_val
        chordVR = vrrotvec(x_hat, chordHats(i,:)) ; 
        rotM_chord_array(:,:,i) = vrrotvec2mat(chordVR) ; 
    end
    [chordHats, ignoreFramesChord] = iterInterpWingVecs(chordHats, ...
        dotThreshChord, rotM_chord_array, N_iter_max, ignoreFrames, ...
        interpType, spanHats) ;
    
    % update data structure
    data_out.([wing_strs{ws} 'ChordHats']) = chordHats ; 
    
    fprintf('Chord correction complete for %s wing \n', wing_strs{ws})
    
    % -------------------------------------------------------
    %% update set of frames to ignore for angle calculations
    if isempty(ignoreFrames_new)
        ignoreFrames_new = sort(unique([ignoreFramesSpan, ...
            ignoreFramesChord]));
    else
        ignoreFrames_new = sort(unique([ignoreFrames_new, ignoreFramesSpan, ...
            ignoreFramesChord]));
    end
end

data_out.ignoreFrames = ignoreFrames_new ;
if  bodyFrameFlag
    data_out = bodyToLabFrame(data_out) ; 
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
%% checks dot product of wing vectors in adjacent frames and interpolates
% perpVecs are used to constrain the estimates of refined vectors, i.e. we
% want to ensure that the chord is always perpendicular to span, so we
% enforce that at each update step by setting perpVecs = spanHats
function [wingVecs, stillBadFrames] = iterInterpWingVecs(wingVecs, ...
    dotThresh, rotM_array, N_iter_max, ignoreFrames, interpType, perpVecs)

if ~exist('interpType','var')
    interpType = 'rot' ; % 'rot' or 'spline'
end
if ~exist('perpVecs','var')
    perpVecs = [] ;
end
% x unit vector is the reference for our rotations
x_hat = [1, 0, 0] ;
N_frames = size(wingVecs,1) ;
searchWindow = 1 ;
maxInterpLength = 6 ;
tol = 1e-8 ;

% get dot product between frame i and (i+1) for all frames (except first)
dotCheck = dot(wingVecs(1:(end-1),:), wingVecs(2:end,:), 2) ;
% also find indices with nan values
nan_idx = any(isnan(wingVecs),2) ;
nan_idx = nan_idx(2:end) ;

% combine these to find indices we'd like to interpolate through
bad_vec_idx = (dotCheck < dotThresh) | nan_idx ;
bad_vec_idx(ignoreFrames) = true ;
% if we have one "good" point sprinkled among bad, probably best to
% interpolate
% for i = (1 + searchWindow) : (length(dotCheck) - searchWindow)
%     i1 = i - searchWindow ;
%     i2 = i + searchWindow ;
%     if (sum(bad_vec_idx(i1:i2)) >= 2*searchWindow)
%         bad_vec_idx(i) = true ;
%     end
% end
% iniate stuff for while loop
cc = 0 ; % counter for iterations
iter_diff = Inf ; % "cost function" difference evaluation

while (sum(bad_vec_idx) > 0) && (cc < N_iter_max) && (iter_diff > tol)
    %tic
    cc = cc + 1 ;
    bad_idx_list = idx_by_thresh(bad_vec_idx) ;
    idx_list_lengths = cellfun(@(y) length(y), bad_idx_list) ;
    bad_idx_list = bad_idx_list(idx_list_lengths <= maxInterpLength) ;
    
    switch interpType
        %------------------------------------------------------------------
        % use matlabs rotational trajectory generation (effectively linear)
        %------------------------------------------------------------------
        case 'rot'
            for j = 1:length(bad_idx_list)
                bad_idx = bad_idx_list{j} + 1 ;  % add 1 because we're essentially doing a "diff"
                
                
                % find start and end points for the trajectory we'll fit
                idx_pre = bad_idx(1) - 1 ;
                idx_post = bad_idx(end) + 1 ;
                
                % make sure we're not going out of array bounds
                if (idx_pre < 1) || (idx_post > N_frames)
                    continue
                end
                
                % get rotation matrices corresponding to start and end
                rotM_pre = rotM_array(:,:,idx_pre) ;
                rotM_post = rotM_array(:,:,(idx_post)) ;
                
                tInterval = [idx_pre, idx_post] ; % frames
                tSamples = idx_pre:idx_post ;
                
                % fit trajectories and store new span vectors
                [rotM_curr,~,~] = rottraj(rotM_pre,rotM_post,tInterval,tSamples) ;
                for k = 2:(size(rotM_curr,3)-1)
                    wingVecs(tSamples(k),:) = (rotM_curr(:,:,k)*x_hat')' ;
                end
                
            end
        %------------------------------------------------------------------
        % use smoothing spline to interpolate 3D trajectory of wingVec tip
        %------------------------------------------------------------------
        case 'spline'
            % fit smoothing interpolants to each coordinate
            frames = (1:N_frames)' ;
            good_ind = find(~bad_vec_idx) + 1 ;
            c_x = fit(frames(good_ind), wingVecs(good_ind,1),...
                'smoothingspline') ; %'smoothingspline'
            c_y = fit(frames(good_ind), wingVecs(good_ind,2),...
                'smoothingspline') ;
            c_z = fit(frames(good_ind), wingVecs(good_ind,3),...
                'smoothingspline') ;
            
            % evaluate splines
            wingVec_interp_x = c_x(frames(2:end)) ;
            wingVec_interp_y = c_y(frames(2:end)) ;
            wingVec_interp_z = c_z(frames(2:end)) ;
            
            wingVecs_interp = [wingVec_interp_x, wingVec_interp_y, ...
                wingVec_interp_z] ;
            wingVecs = [wingVecs(1,:) ; wingVecs_interp] ;
            % normalize
            wingVecs = wingVecs ./ repmat(myNorm(wingVecs),1,3) ;
        otherwise
            keyboard
    end
    
    % --------------------------------------
    % enforce orthogonality, if applicable
    if ~isempty(perpVecs)
        wingVecsPerp = wingVecs - perpVecs.*(dot(wingVecs,perpVecs,2)) ;
        wingVecsPerp = wingVecsPerp./myNorm(wingVecsPerp) ;
        wingVecs = wingVecsPerp ;
    end
    
    % ----------------------------------------------------------------
    % check adjacent frame dot product for newly calculated wing vecs
    dotCheck_new = dot(wingVecs(1:(end-1),:), wingVecs(2:end,:), 2) ;
    iter_diff = nanmean(1-dotCheck) - nanmean(1-dotCheck_new) ;
    dotCheck = dotCheck_new ;
    bad_vec_idx = (dotCheck < dotThresh) ;
    % fill holes in bad index list
    for i = (1 + searchWindow) : (length(dotCheck) - searchWindow)
        i1 = i - searchWindow ;
        i2 = i + searchWindow ;
        if (sum(bad_vec_idx(i1:i2)) >= 2*searchWindow)
            bad_vec_idx(i) = true ;
        end
    end
    %toc
    %disp(sum(1-dotCheck))
end
disp(cc)
stillBadFrames = intersect((find(bad_vec_idx) + 1)', ignoreFrames) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SNIPPETS THAT MAY BE USEFUL LATER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     bad_span_idx(1:searchWindow) = false ; 
%     bad_span_idx((end-searchWindow):end) = false ; 
%     