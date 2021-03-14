% -------------------------------------------------------------------------
% function to correct errors in span (and chord?) estimation by looking at
% the time series of vector orientations. May need to use images/voxels as
% well
% -------------------------------------------------------------------------
function data_out = refineWingVecIter_v2(data_in, interpType, bodyFrameFlag)
% ------------------------
%% params and inputs
if ~exist('interpType','var') || isempty(interpType)
    interpType = 'spline' ; % 'spline' or 'rot'
end
if ~exist('bodyFrameFlag','var') || isempty( bodyFrameFlag)
    bodyFrameFlag = false ; % check wing vectors in body frame of reference
end

% minimum value for dot product wing vector combinations
dotThreshSpan = 0.90 ;
dotThreshChord = 0.2 ; 
dotThreshChordFlip = 0.85 ;
dotThreshChordStrict = 0.999 ;
dotThreshChordVel = 0.20 ;
% minimum number of indices between estimated flip times
minFlipInterval = 12 ;
flipPad = 3 ; % for a given flip time t, we say (t-flipPad):(t+flipPad) is when the wing is going "slow"
% method of interpolation
% max number of iterations to converge on good set of wing vecs
N_iter_max = 20 ;
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
    % get identifiers for the wing we're dealing with
    wing_str = wing_strs{ws} ; 
    wing_str_cap = regexprep(lower(wing_str),'(\<[a-z])','${upper($1)}') ;
    wing_sum_str = wing_num_strs{ws} ; 
    
    % read in wing vectors from data structure
    spanHats = data_out.([wing_str 'SpanHats']) ;
    chordHats = data_out.([wing_str 'ChordHats']) ;
    chordAltHats = data_out.(['chord' wing_sum_str 'AltHats']) ;
    %... also get diagonals (for evaluating chord)
    diag1 = data_out.(['diag' wing_num_strs{ws} '1' wing_str_cap]) ; 
    diag2 = data_out.(['diag' wing_num_strs{ws} '2' wing_str_cap]) ; 
    
    % seems like chordAltHats are not always perpendicular to span...
    chordAltHats = chordAltHats - ...
        repmat(dot(chordAltHats,spanHats,2),1,3).*spanHats ;
    chordAltHats = chordAltHats./repmat(myNorm(chordAltHats),1,3) ;
    
    % depending on whether we're on the left or right side, the midstroke
    % point (90 deg) is the +/- y axis (left/right, respectively)
    if strcmp(wing_strs{ws},'right')
        spanFrameRef = [0, -1, 0] ;
    elseif strcmp(wing_strs{ws},'left')
        spanFrameRef = [0, 1, 0] ;
    else
        disp('wing side can only be right or left')
    end
    
    % ---------------------------------------------------------------------
    %% initialize new data storage 
    newChordHats = chordHats ; 
    newChordAltHats = chordAltHats ; 
    newDiag1 = diag1 ; 
    newDiag2 = diag2 ; 
    
    % in case we're correction the span using matlab's function to
    % interpolate rotation matrices (or quaternions), get rotation matrices
    % that correspond to span positions
    if strcmp(interpType, 'rot')
        rotM_span_array = zeros(3,3,N_frames) ;
        span_nan_idx = any(isnan(spanHats(2:end,:)),2) ;
        loop_val = find(~span_nan_idx)' + 1 ;
        for i = loop_val
            spanVR = vrrotvec(x_hat, spanHats(i,:)) ;
            rotM_span_array(:,:,i) = vrrotvec2mat(spanVR) ;
        end
    else
        rotM_span_array = [] ; 
    end

    % ---------------------------------------------
    %% first deal with spans
    fprintf('Correcting spans for %s wing... \n', wing_strs{ws})
    % find frames where dot product between adjacent frames is off and
    % interpolate through them. do this until all frames look good or we
    % reach iteration limit
    [newSpanHats, ignoreFramesSpan, spanVel] = iterInterpWingVecs(spanHats, ...
        dotThreshSpan, rotM_span_array, N_iter_max, ignoreFrames,interpType) ;
    
    % add new spans to data structure
    data_out.([wing_strs{ws} 'SpanHats']) = newSpanHats ;
    
    fprintf('Span correction complete for %s wing \n', wing_strs{ws})
    
    % --------------------------------------------------------
    %% next check chords
    % going to separate frames into "fast" and "slow" frames. In "fast"
    % frames, chord should align with wing velocity. in slow frames. we're
    % just looking for a rotation with a consistent sign
    fprintf('Correcting chords for %s wing... \n', wing_strs{ws})
    
    % first make sure that chords are perpendicular to new span
    chordHats = chordHats - dot(chordHats,newSpanHats,2).*newSpanHats ;
    chordHats = chordHats./myNorm(chordHats) ;
    chordAltHats = chordAltHats - ...
        dot(chordAltHats,newSpanHats,2).*newSpanHats ;
    chordAltHats = chordAltHats./myNorm(chordAltHats) ;
    
    % ------------------------
    %% deal with "fast" case
    
    swap_flag = false(N_frames, 1) ; 
    % compare the direction of the chord against the direction of the
    % (hopefully now corrected) motion of the span vector
%     spanHatVel_xy = [spanHatVel(:,1:2), zeros(N_frames,1)] ;
%     spanHatVel_xy = spanHatVel_xy./myNorm(spanHatVel_xy) ;
    spanVelHat = spanVel./myNorm(spanVel) ; 
    spanVel_norm = myNorm(spanVel) ;
    
    % get chord and alt chord direction in xy plane
%     chordHats_xy = [chordHats(:,1:2), zeros(N_frames,1)] ;
%     chordHats_xy = chordHats_xy./myNorm(chordHats_xy) ;
%     
%     chordAltHats_xy = [chordAltHats(:,1:2), zeros(N_frames,1)] ;
%     chordAltHats_xy = chordAltHats_xy./myNorm(chordAltHats_xy) ;
    
    % stroke plane normal vectors--for now just taking z axis, but could
    % improve in future (affine fit to chord vectors?)
    strokePlaneNormals = repmat([0,0,1],N_frames,1) ; 
    
    % compile chordHat options and dot products with xy span velocity
%     chordHatOptions = cat(3,chordHats, chordAltHats) ;
%     chordDot_xy = [dot(spanHatVel_xy, chordHats_xy,2), ...
%         dot(spanHatVel_xy, chordAltHats_xy,2)] ;
%     % also get "dot product" with positive z direction
%     chordDot_z = [chordHats(:,3), chordAltHats(:,3)] ; 
    
    % ---------------------------------------------------------------------
    % get a rough estimation of flip times. use this to  identify periods
    % where the wing is moving "fast"--in these cases, chord should align
    % with wing velocity
    [~, estFlipTimes] = findpeaks(-1*spanVel_norm,...
        'MinPeakDistance',minFlipInterval) ;
    flip_idx = false(N_frames,1) ;
    flip_idx(estFlipTimes) = true ;
    flip_idx = bwmorph(flip_idx,'thicken',flipPad) ;
    flip_idx = circshift(flip_idx, flipPad-round(flipPad/2)) ;
    
    % for the non-flipping times, take the chord that matches the span
    % velocity best (also check z component
    % NEED TO COMPARE DIFFERENCE--IF THEY'RE BOTH ROUGHLY IN THE RIGHT
    % DIRECTION, NEED OTHER CHECK
    %[~ , chord_sort_idx] = sort(chordDot_xy,2, 'descend') ;
    %[chordDot_xy_sort, sort_idx] = sort(chordDot_xy,2,'descend') ;
    fast_ind = find(~flip_idx) ;
    for it = fast_ind'
        c1 = chordHats(it, :) ;
        c2 = chordAltHats(it,:) ;
        %c1_xy = chordHats_xy(it, :) ;
        %c2_xy = chordAltHats_xy(it, :) ;
        
        v = spanVelHat(it,:) ;
        
        % if needed, invert c1 and c2 such that they point in the direction of
        % tip velocity
        
        dot_vel1 = dot(c1, v) ; 
        if (dot_vel1 < 0)
            c1 = - c1 ;
            dot_vel1 = dot(c1, v) ; 
        end
        dot_vel2 = dot(c2, v) ; 
        if (dot_vel2 < 0)
            c2 = - c2 ;
            dot_vel2 = dot(c2, v) ; 
        end
        
        % if there is only one chord vector that points up w.r.t. stroke plane,
        % choose this one
        dot1 = dot(c1, strokePlaneNormals(it,:)) ;
        dot2 = dot(c2, strokePlaneNormals(it,:)) ;
        
        % make various checks to see if we should swap
        strokePlaneCheck = (dot2>0) && (dot1<0) ; 
        diagCheck = (dot2>0) && (diag2(it)/diag1(it)>1.2) ; 
        %velCheck = (dot2>0) && ((dot_vel2 - dot_vel1) > dotThreshChordVel) ; 
        if (strokePlaneCheck) || (diagCheck) %|| (velCheck)
            swap_flag(it)     = true ;
            newChordHats(it,:)    = c2 ; % swap
            newChordAltHats(it,:) = c1 ;
            newDiag1(it)      = diag2(it) ;
            newDiag2(it)      = diag1(it) ;
        else
            newChordHats(it,:)    = c1 ; % do not swap
            newChordAltHats(it,:) = c2 ;
            newDiag1(it)      = diag1(it) ;
            newDiag2(it)      = diag2(it) ;
        end
        %         chordHatOptions(f_ind,:,:) = ...
        %             chordHatOptions(f_ind, :, chord_sort_idx(f_ind,:)) ;
        %         chordHats(f_ind,:) = ...
        %             squeeze(chordHatOptions(f_ind,:,chord_sort_idx(f_ind, 1))) ;
        %         chordAltHats(f_ind,:) = ...
        %             squeeze(chordHatOptions(f_ind,:,chord_sort_idx(f_ind, 2))) ;
        
        %     chordHats = squeeze(chordHatOptions(:,:,1)) ;
        %     chordAltHats = squeeze(chordHatOptions(:,:,2)) ;
    end
    
    % find the mean diag value of the frames processed so far
    ind = (newDiag1>0) ; % values of 0 correspond to cases where wing was not identified...
    meanDiag = mean(newDiag1(ind)) ;
    stdDiag  = std(newDiag1(ind)) ;
    
    % check if any of these newly-assgned chords have negative z component
    zCheck = (newChordHats(:,3) < 0) ;
    newChordHats(zCheck & ~ flip_idx,:) = -1*newChordHats(zCheck & ~flip_idx,:) ;
    zCheckAlt = (newChordAltHats(:,3) < 0) ;
    newChordAltHats(zCheckAlt & ~ flip_idx,:) = ...
        -1*newChordAltHats(zCheckAlt & ~flip_idx,:) ;
    
    % find any bad "fast" frames and try to interpolate through those
%     chordHatsTmp_xy = [newChordHats(:,1:2), zeros(N_frames,1)] ;
%     chordHatsTmp_xy = chordHatsTmp_xy ./myNorm(chordHatsTmp_xy) ;
    chordVelDot = dot(newChordHats , spanVelHat,2) ;
    
    bad_fast_idx = (chordVelDot < dotThreshChordVel) & (~flip_idx) ;
    nan_idx = any(isnan(newChordHats),2) ;
    %bad_fast_frames = find(bad_fast_idx) ;
    chordHatsInterp = splineWingHatVec(newChordHats, ~(bad_fast_idx | nan_idx),...
        'cubicinterp') ; % 'cubicinterp'
    newChordHats(bad_fast_idx,:) = chordHatsInterp(bad_fast_idx,:) ;
    
    % ---------------------------------------------
    %% now deal with "slow" frames
    % get matrices that rotate span vector onto yhat (or -1*yhat, if right
    % side). in this reference frame, chords should just rotate in xz plane
    rotM_spanFrame = nan(3,3,N_frames) ;
    chordHatSpanFrame = nan(size(newChordHats)) ;
    chordAltHatSpanFrame = nan(size(newChordAltHats)) ;
    spanRotCheck = nan(size(newSpanHats)) ;
    phiList = nan(N_frames,1) ;
    for ind = 1:N_frames
        % get rotation matrix by first undoing yaw and pitch (stroke and
        % deviation, in wing nomenclature) then rotating from x axis to +y
        % (or -y) axis
        phi = atan2(newSpanHats(ind,2), newSpanHats(ind,1))  ;
        theta = asin(newSpanHats(ind,3));
        rotM1 = eulerRotationMatrix(phi, theta, 0) ;
        rotM2 = eulerRotationMatrix(sign(spanFrameRef(2))*(-1*pi/2),0,0) ;
        rotM = rotM2*rotM1 ;
        
        % apply rotations to chord, alt chord, and span. also store stroke
        % angle to check flip direction
        chordHatSpanFrame(ind,:) = (rotM*newChordHats(ind,:)')' ;
        chordAltHatSpanFrame(ind,:) = (rotM*newChordAltHats(ind,:)')' ;
        spanRotCheck(ind,:) = (rotM*newSpanHats(ind,:)')' ;
        rotM_spanFrame(:,:,ind) = rotM ;
        phiList(ind) = phi ;
    end
    
    % check for magnitude and direction of rotation between subsequent
    % frames (in window around slow frames)
    slow_bouts = idx_by_thresh(flip_idx) ;
    bad_slow_idx = false(N_frames, 1) ;
    for indd = 1:length(slow_bouts)
        bout = slow_bouts{indd} ;
        i1 = max([1, bout(1) - 3]) ;
        i2 = min([N_frames, bout(end) + 3]) ;
        
        if (i1 == 1)
            % haven't figured out what to do when we have a bad estimate of
            % chord in first frame yet
            continue
        end
        
        % check expected sign of rotation using stroke angle
        meanPhi = nanmean(phiList(bout)) ;
        if (meanPhi*sign(spanFrameRef(2)) < pi/2)
            rotSignPhi = -1 ;
        elseif (meanPhi*sign(spanFrameRef(2)) > pi/2)
            rotSignPhi = 1 ;
        else
            disp('bad estimate of phi?')
            keyboard
        end
        % check full angle sweep from i1 to i2
        dotFull = dot(chordHatSpanFrame(i1,:), ...
            chordHatSpanFrame(i2,:));
        rotAngAvg = acos(dotFull)/length(bout) ; 
        crossFull = cross(chordHatSpanFrame(i1,:), ...
            chordHatSpanFrame(i2,:)) ;
        rotSign = sign(crossFull(2)) ;
        
        % loop through and check framewise rotation
        for i = (i1+2):(i2-3)
            % get dot and cross products
            dotCurr = dot(chordHatSpanFrame(i,:),chordHatSpanFrame(i+1,:)) ;
            dotAlt = dot(chordHatSpanFrame(i,:), ...
                chordAltHatSpanFrame(i+1,:)) ;
            %disp([num2str(dotCurr) ' ' num2str(dotAlt)])
            crossCurr = cross(chordHatSpanFrame(i,:), ...
                chordHatSpanFrame(i+1,:)) ;
            crossAlt = cross(chordHatSpanFrame(i,:), ...
                chordAltHatSpanFrame(i+1,:)) ;
            
            % make checks to see if the either the original or alternate
            % chord have rotations that make sense. if neither, flag for
            % interpolation
            checkCurr = ((sign(crossCurr(2)) == rotSign) | ...
                (dotCurr > dotThreshChordStrict)) & (dotCurr >= dotThreshChordFlip) ;
            checkAlt = ((sign(crossAlt(2)) == rotSign) | ...
                (dotAlt > dotThreshChordStrict)) & (dotAlt >= dotThreshChordFlip) ;
            % if the original chord gives a rotation in the correct
            % direction or just doesn't change much from the previous
            % frame, keep it. if that's not the case, but altChord fits
            % that, take it instead. otherwise, flag that frame
            % * NB: we're looking at the rotation from i to i+1, and maybe
            % changing the chord in frame i+1
            if ~checkCurr && checkAlt
                % swap chords
                tmp = chordHatSpanFrame(i+1,:) ;
                chordHatSpanFrame(i+1,:) = chordAltHatSpanFrame(i+1,:) ;
                chordAltHatSpanFrame(i+1,:) = tmp ;
            elseif ~checkCurr && ~checkAlt
                bad_slow_idx(i+1) = true ;
                % also, to not mess up subsequent checks, make a rough
                % guess of where this chord 'should' be based on average
                % angular velocity of flip
%                 rot_ang_vel = acos(dot(chordHatSpanFrame(i-1,:),...
%                     chordHatSpanFrame(i,:))) ;
                temp_rotM = axang2rotm([spanFrameRef, ...
                    sign(spanFrameRef(2))*rotAngAvg*rotSign]) ;
                chordHatSpanFrame(i+1,:) = ...
                    (temp_rotM*chordHatSpanFrame(i,:)')' ;
                
            end
        end
       
        % ------------------------------------------
        % plot to visualize chords during flip
        if (0)
            figure ;
            hold on
            N_vecs = length(i1:i2) ;
            ind_start = i1 ;
            cmat = brewermap(N_vecs, 'RdYlBu') ;
            cc = 1 ;
            for j = (ind_start):(ind_start+N_vecs-1)
                plot([0, chordHatSpanFrame(j,1)], [0, chordHatSpanFrame(j,3)], ...
                    'Color',cmat(cc,:),'LineWidth',2.5)
                %             plot([0, chordAltHatSpanFrame(j,1)], [0, chordAltHatSpanFrame(j,3)], ...
                %                 'Color',cmat(cc,:),'LineWidth',2.5,'LineStyle',':')
                text(chordHatSpanFrame(j,1) - 0.01,...
                    chordHatSpanFrame(j,3) + 0.03,...
                    num2str(j))
                %             text(chordAltHatSpanFrame(j,1) + 0.01, ...
                %                 chordAltHatSpanFrame(j,3) + 0.03,...
                %                 [num2str(j) '^~'])
                cc = cc + 1 ;
            end
            %axis tight
            axis equal
            grid on
        end
        
    end
    
    % -------------------------------
    % interpolate bad frames
    nan_idx = any(isnan(chordHatSpanFrame),2) ;
    chordHatSpanFrame = splineWingHatVec(chordHatSpanFrame, ...
        ~(bad_slow_idx | nan_idx), 'cubicinterp') ; 
%     chordAltHatSpanFrame = splineWingHatVec(chordAltHatSpanFrame, ~bad_slow_idx, ...
%         'cubicinterp') ; 
    
    
    % ---------------------------------------
    % transform chords back into body frame
    for ind = 1:N_frames
        rotM = squeeze(rotM_spanFrame(:,:,ind)) ;
        newChordHats(ind,:) = (rotM'*chordHatSpanFrame(ind,:)')' ;
        newChordAltHats(ind,:) = (rotM'*chordAltHatSpanFrame(ind,:)')' ;
    end
    
    % make sure that chords are perpendicular to spans
    newChordHats = newChordHats - ...
        dot(newChordHats,newSpanHats,2).*newSpanHats ;
    newChordHats = newChordHats./myNorm(newChordHats) ;
    newChordAltHats = newChordAltHats - ...
        dot(newChordAltHats,newSpanHats,2).*newSpanHats ;
    newChordAltHats = newChordAltHats./myNorm(newChordAltHats) ;
    
    %update data structure
%     data_out.([wing_strs{ws} 'ChordHats']) = newChordHats ;
%     data_out.(['chord' wing_num_strs{ws} 'AltHats']) = newChordAltHats ;
%     
%     %disp(sum(1-chordDot))
    ignoreFramesChord = [] ;
    %-------------------------------------------------------------------
    %% then use same trick with span to interpolate through bad points
    % get (chord) matrices for each frame
    %     fprintf('Correcting chords for %s wing... \n', wing_strs{ws})
%     if strcmp(interpType, 'rot')
%         rotM_chord_array = zeros(3,3,N_frames) ;
%         chord_nan_idx = any(isnan(chordHats(2:end,:)),2) ;
%         loop_val = find(~chord_nan_idx)' + 1 ;
%         for i = loop_val
%             chordVR = vrrotvec(x_hat, chordHats(i,:)) ;
%             rotM_chord_array(:,:,i) = vrrotvec2mat(chordVR) ;
%         end
%     else
%         rotM_chord_array = [] ;
%     end
%     
%     [newChordHats, ignoreFramesChord] = iterInterpWingVecs(newChordHats, ...
%         dotThreshChord, rotM_chord_array, 3, ignoreFrames, ...
%         interpType, newSpanHats) ;
    
    % update data structure
    data_out.([wing_strs{ws} 'ChordHats']) = newChordHats ;
    data_out.(['chord' wing_num_strs{ws} 'AltHats']) = newChordAltHats ;
    
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
function [wingVecs, stillBadFrames, wingVecVel] = iterInterpWingVecs(wingVecs, ...
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
tol = 5e-4 ;

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
            %frames = (1:N_frames)' ;
            good_ind = find(~bad_vec_idx) + 1 ;
            [wingVecs, wingVecVel,~] = splineWingHatVec(wingVecs, good_ind, ...
                'smoothingspline') ;
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
    
    % --------------------------------------------
    % get velocities associated with wing vec
    if ~strcmp(interpType,'spline')
        disp('under construction')
        wingVecVel = [] ;
    end
    % ----------------------------------------------------------------
    % check adjacent frame dot product for newly calculated wing vecs
    dotCheck_new = dot(wingVecs(1:(end-1),:), wingVecs(2:end,:), 2) ;
    iter_diff = nanmean(1-dotCheck) - nanmean(1-dotCheck_new) ;
    %disp(iter_diff)
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

% -------------------------------------------------------------------------
%% use a spline to interpolate wing vectors
function [wingVecSpline, wingVecVel, wingVecAccel] = ...
    splineWingHatVec(wingVecs, interp_ind, splineType)

% get frames
frames = (1:size(wingVecs,1))' ;

% fit splines
c_x = fit(frames(interp_ind), wingVecs(interp_ind,1), splineType) ; %'smoothingspline'
c_y = fit(frames(interp_ind), wingVecs(interp_ind,2), splineType) ;
c_z = fit(frames(interp_ind), wingVecs(interp_ind,3), splineType) ;

% evaluate splines
wingVec_interp_x = c_x(frames(2:end)) ;
wingVec_interp_y = c_y(frames(2:end)) ;
wingVec_interp_z = c_z(frames(2:end)) ;

% get smoothed result
wingVecSpline = [wingVec_interp_x, wingVec_interp_y, ...
    wingVec_interp_z] ;
wingVecSpline = [wingVecs(1,:) ; wingVecSpline] ;
% normalize
wingVecSpline = wingVecSpline ./ repmat(myNorm(wingVecSpline),1,3) ;

if nargout > 1
    [vx, ax] = differentiate(c_x, frames) ;
    [vy, ay] = differentiate(c_y, frames) ;
    [vz, az] = differentiate(c_z, frames) ;
    
    wingVecVel = [vx, vy, vz] ;
    wingVecAccel = [ax, ay, az] ;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SNIPPETS THAT MAY BE USEFUL LATER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     bad_span_idx(1:searchWindow) = false ;
%     bad_span_idx((end-searchWindow):end) = false ;
%