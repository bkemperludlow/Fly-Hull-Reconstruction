% -------------------------------------------------------------------------
% my attempt at a function to fix issues with chord vector, based on
% Tsevi's "refineChordVector.m", but, since his wasn't completed and I
% don't fully understand what he was driving at, I'm going to have to
% re-work it so it makes sense to me
%
%{

% Tsevi's notes for his function:

    -find the frames that correspond to forward and back stroke
    -find the stroke plane
    -still need to handle cases where the wing points down with respect to
     the stroke plane
    -also need to handle cases where "diag" difference is very large.
%}
% MY TO-DO:
%   - check span for errors/interpolate through
% -------------------------------------------------------------------------
function data_out = refineChordVectorSam(data, wingSide, ...
    largePertFlag, SPEED_THRESHOLD, DELTA, plotFlag)
% -------------------------------------------------------------------
%% inputs and params
if ~exist('largePertFlag', 'var') || isempty(largePertFlag)
    largePertFlag = false ; % is this a large perturbation?
end
if ~exist('SPEED_THRESHOLD', 'var') || isempty(SPEED_THRESHOLD)
    SPEED_THRESHOLD = 1 ; % m/s
end
if ~exist('DELTA', 'var') || isempty(DELTA)
    DELTA = 17; % take +17 and -17 points from each time point. to find the "current" stroke plane.
end
if ~exist('plotFlag', 'var') || isempty(plotFlag)
    plotFlag = false ;
end

debugFlag1 = false ; % check out wing voxels and chord in case of bad chord in fast pass 1
smoothSpanFlag = false ; % try to smooth/interp span at outset?
% ---------------------------------------
% make sure wing side is in proper format
if strcmpi(wingSide, 'right') || strcmpi(wingSide,'R')
    wingSide = 'right' ;
    wingSideNum = '1' ;
    otherSide = 'left' ;
elseif strcmpi(wingSide, 'left') || strcmpi(wingSide,'L')
    wingSide = 'left' ;
    wingSideNum = '2' ;
    otherSide = 'right' ;
else
    fprintf('Invalid wingSide selection: %s \n', wingSide)
    keyboard
end

% -------------------------------------
% moving slope derivative params
movslope_order = 2 ;
movslope_len = max([8, movslope_order+1]) ; % make sure it's larger than order

% -------------------------------------------------------------------------
%% read data from structure
% first transform data to body frame:
data_bodyFrame = labToBodyFrame(data, largePertFlag) ;

% get interpolated wing tip and wing cm trajectory
[wingCM, ~, wingTip] = interpolateWingCM(data_bodyFrame, wingSide) ;
[~, ~, wingTipOther] = interpolateWingCM(data_bodyFrame, otherSide) ;

% convert wing tip/cm to meters
voxelSize = data_bodyFrame.params.voxelSize ;
wingTip      = voxelSize * wingTip ;
wingCM       = voxelSize * wingCM ;
wingTipOther = voxelSize * wingTipOther ;

% also read out wing vectors:
span = data_bodyFrame.([wingSide 'SpanHats']) ;
mainChord = data_bodyFrame.([wingSide 'ChordHats']) ;
altChord = data_bodyFrame.(['chord' wingSideNum 'AltHats']) ;
diag1 = data_bodyFrame.(['diag' wingSideNum '1' upper(wingSide(1)) ...
    wingSide(2:end)]) ;
diag2 = data_bodyFrame.(['diag' wingSideNum '2' upper(wingSide(1)) ...
    wingSide(2:end)]) ;

% body axis vector
AHat = data_bodyFrame.AHat ;

% other params
dt      = 1 / data_bodyFrame.params.fps ;
N       = size(wingTip,1) ;
tvec    = (0:N-1)' * dt ;
tvec_ms = tvec * 1000 ;

% -----------------------------------------------------------------------
%% initialize outputs
newChord     = zeros(N, 3) ;
newAltChord  = zeros(N, 3) ;
newDiag1     = zeros(N, 1) - 1; % a negative value indicates a frame we did not get to.
newDiag2     = zeros(N, 1) - 1;

swap_flag = false(N,1) ;
error_flag = false(N,1) ; 

% -----------------------------------------------------------------------
%% smooth/interpolate span?
if smoothSpanFlag
    %span_smooth = nan(size(span)) ;
    for dim = 1:3
        [~, hampel_idx] = hampel(span(:,dim)) ;
        span(hampel_idx,dim) = nan ;
        span(:,dim) = smooth(tvec, span(:,dim), 17, 'rloess') ;
    end
    
    % normalize smoothed span
    span = span./myNorm(span) ;
    if (0)
        figure ;
        for dim = 1:3
            subplot(3,1,dim)
            hold on
            plot(tvec, span(:,dim),'o-')
            %plot(tvec, span_smooth(:,dim),'-')
        end
    end
end
% -----------------------------------------------------------------------
%% find a fit to stroke plane for each wingbeat (based on fixed window)
% this fit should go through wing tip on both sides

strokePlaneNormals = zeros(N,3) ;
for it=1:N
    t1 = max([1, it-DELTA]) ;
    if (t1==1)
        t2 = 2*DELTA+1 ;
    else
        t2 = min([N, it+DELTA]);
        if (t2==N)
            t1 = N-2*DELTA ;
        end
    end
    if (t2-t1~=2*DELTA)
        disp('problem') ;
        keyboard ;
    end
    
    [n, V, p] = affine_fit([wingTip(t1:t2,:) ; wingTipOther(t1:t2,:)]) ;
    strokePlaneNormals(it,:) = n' ;
    
end

% ------------------------------------------------------------------------
%% get tip velocity (in body frame)

% make sure moving window is less than array size
movslope_len = min([movslope_len, N-1]) ;

% loop through x, y, z and get derivatives
tipVel = zeros(size(wingTip)) ;
for dim = 1:3
    tipVel(:,dim) = (1/dt)*movingslope(wingTip(:,dim),...
        movslope_len, movslope_order) ;
end

% calc velocity in stroke plane
tipVel_plane = tipVel - strokePlaneNormals.*dot(tipVel, strokePlaneNormals, 2)  ;

% also get projection of body axis in stroke plane
AHat_plane = AHat - strokePlaneNormals.*dot(AHat, strokePlaneNormals, 2)  ;
AHat_plane = AHat_plane./myNorm(AHat_plane) ;

% also get wing tip total speed
speed = myNorm(tipVel) ;
speed_plane = myNorm(tipVel_plane) ;

% calc component of tip velocity along the ahat, both vectors in the stroke plane
vv       = dot(tipVel_plane, AHat_plane,2) ;
fast_ind = ( abs(vv) > SPEED_THRESHOLD) ; % determine "fast" criterion
% ------------------------------------------------------------------------
%% deal with both "fast" and "slow" slow frames (differently)
% if the wing is moving "fast", then choose the chord that is more aligned
% with the wing-tip velocity

% if the wing is moving "slow" - this is the harder case. dothe following:
%
%{
find a time range around a wing flip
guess when the flip was, i.e. in the middle between two adjacent frames. if
frames are f and f+1, then assume flip was between them

before the flip chord vector should point consistently in one direction
after the flip the chord vector should point consistently in the opposite
direction. "direction" here is with respect to the tip velocity vector in
the stroke plane.

the "before" direction is the wing-tip velocity in the beginning of the
time interval

the "after" direction is the wing-tip velocity in the end of the time
interval
%}

% -------------------------------------------------------
% First handle the frames with higher wing-tip velocity
% --> find chord more aligned with wing tip velocity
for it= 1:N
    if (~fast_ind(it))
        continue ;
    end
    v = tipVel_plane(it,:) ; % tip_v(it,:) ; % curret tip velocity
    v = v / norm(v) ; % normalize to get direction only
    
    c1 = mainChord(it,:) ;
    c2 = altChord(it,:) ;
    
    % if needed, invert c1 and c2 such that they point in the direction of
    % tip velocity
    
    if (dot(c1, v)<0)
        c1 = - c1 ;
    end
    
    % if the main chord doesn't point up w.r.t. stroke plane, recalculate
    dot1 = dot(c1, strokePlaneNormals(it,:)) ;
    
    if ( dot1 < 0 )
        if debugFlag1
            h_chord = showChordAndVox(data_bodyFrame, wingSide, it, c1, ...
                c2, span(it,:), (1/voxelSize)*wingTip(it,:), ...
                tipVel_plane(it,:))  ;
            keyboard
        end
        
        [chord_new, altChord_new, diag1_new, diag2_new, error_flag(it)] = ...
            reCalcChord(data_bodyFrame, wingSide, it, c1, c2, ...
            span(it,:), tipVel_plane(it,:), strokePlaneNormals(it,:),...
            'fast')  ;
        
        swap_flag(it)     = true ;
        newChord(it,:)    = chord_new ; % swap
        newAltChord(it,:) = altChord_new ;
        newDiag1(it)      = diag1_new ;
        newDiag2(it)      = diag2_new ;
    else
        newChord(it,:)    = c1 ; % do not swap
        newAltChord(it,:) = c2 ;
        newDiag1(it)      = diag1(it) ;
        newDiag2(it)      = diag2(it) ;
    end
end

% find the mean diag value of the frames processed so far
ind = (newDiag1>0) ; % values of 0 correspond to cases where wing was not identified...
meanDiag = mean(newDiag1(ind)) ;
stdDiag  = std(newDiag1(ind)) ;

% ----------------------------------------------------------------------
%% pass 2
% go over all fast frames again, and check for the ones with diag1 too
% small and diag2 around the "good size"
for it=1:N
    if (~fast_ind(it))
        continue ;
    end
    if (newDiag1(it) < 0.55*newDiag2(it)) &&  ... % if diag1 is smaller than half diag2
            (abs(newDiag2(it)-meanDiag)/stdDiag < 3) % and if diag2 is with 3 sigma from mean diag
        % then swap
        swap_flag(it) = ~swap_flag(it) ;
        
        tmp = newChord(it,:) ;
        newChord(it,:) = newAltChord(it,:) ;
        newAltChord(it,:) = tmp;
        
        tmp =  newDiag1(it) ;
        newDiag1(it) = newDiag2(it) ;
        newDiag2(it) = tmp ;
        %disp(it);
        clear tmp
    end
end

% check the histogram again
ind = (newDiag1>0) ; % values of 0 correspond to cases where wing was not identified...
meanDiag = mean(newDiag1(ind)) ;
stdDiag  = std(newDiag1(ind)) ;

% -------------------------------------------------------------------------
%% pass 3
% look at chunks of fast moving frames and make sure chord is moving in a
% reasonable way

% get chunks of "fast" frames
fast_bouts = idx_by_thresh(fast_ind) ;

% also initialize list of body to (unpitched) wing frame rotation matrices
body2WingArray = nan(N, 3, 3) ;
newChordRot = nan(size(newChord)) ;
newAltChordRot = nan(size(newAltChord)) ;

% loop through fast bouts
for bout_ind = 1:length(fast_bouts)
    frame_ind = fast_bouts{bout_ind} ;
    
    % -------------------------------------------------------
    % first loop -- get rotation matrices and rotate chord
    for it = frame_ind'
        % get body to (unpitched) wing frame rotation matrix
        phi_w = atan2(span(it, 2), span(it, 1)) ; % wing stroke angle
        theta_w = asin(span(it, 3)) ; % wing deviation angle
        rotM = body2WingRotMat(phi_w, theta_w, 0) ; % rotation mat
        
        % rotate chords and store rotation mat
        newChordRot(it,:) = (rotM*newChord(it,:)')' ;
        newAltChordRot(it,:) = (rotM*newAltChord(it,:)')' ;
        body2WingArray(it,:,:) = rotM ;
    end
    
    % plot results?
    if (0)
        figure ;
        hold on
        colorMat = parula(length(frame_ind)) ;
        cc = 1 ;
        for itt = frame_ind'
            plot([0, newChordRot(itt, 2)], [0, newChordRot(itt, 3)], 'o-', ...
                'Color', colorMat(cc,:),'LineWidth',1.25)
            plot([0, newAltChordRot(itt, 2)], [0, newAltChordRot(itt, 3)], ':', ...
                'Color', colorMat(cc,:),'LineWidth',0.75)
            cc = cc + 1 ;
        end
        axis equal
        box on
        grid on
        
    end
    
    % ----------------------------------------------
    % check for small number of obvious outliers
    chordRotMat = newChordRot(frame_ind,:) ;
    cosDistMat = squareform(pdist(chordRotMat,'cosine')) ;
    medianCosDist = median(cosDistMat,2) ;
    % find chords that tend to be different from many others
    bad_idx = (medianCosDist > 0.3) ;
    bad_ind = find(bad_idx) ;
    % check that this is only a few--if so, first try swapping them if that
    % doesn't work
    if (numel(bad_ind) < (length(frame_ind)/3))
        for ind = bad_ind'
            indd = frame_ind(bad_ind) ;
            altChordRotCurr = newAltChordRot(indd,:) ;
            distCheck = median(pdist2(chordRotMat(~bad_idx,:), ...
                altChordRotCurr, 'cosine')) ;
            % if the alt chord fits much better, add this to bout matrix and
            % perform swap
            if (distCheck < 0.3)
                % add to bout matrix
                chordRotMat(bad_ind,:) = altChordRotCurr ;
                
                % -----------------------------
                % do swap
                swap_flag(indd) = ~swap_flag(indd) ;
                
                % swap vectors
                tmp = newChord(indd,:) ;
                newChord(indd,:) = newAltChord(indd,:) ;
                newAltChord(indd,:) = tmp ;
                
                % swap rotated vectors
                tmp2 = newChordRot(indd,:) ;
                newChordRot(indd,:) = newAltChordRot(indd,:) ;
                newAltChordRot(indd,:) = tmp2 ;
                
                % swap diagonals
                tmp3 = newDiag1(indd) ;
                newDiag1(indd) = newDiag2(indd) ;
                newDiag2(indd) = tmp3 ;
            else
                % if we can't solve this with a swap, we'll have to
                % interpolate, so for now just make this value nan (we'll do
                % interpolation in next pass)
                chordRotMat(bad_ind,:) = nan ;
                newChord(indd,:) = nan ;
                newDiag1(indd) = nan ;
                error_flag(it) = true ;
            end
        end
    end
    
    % -----------------------------------------
    % interpolate through any nan points
    nan_idx = any(isnan(chordRotMat),2) ;
    if (sum(~nan_idx) < 3)
        continue
    elseif (sum(nan_idx) > 0)
        x = 1:size(chordRotMat,1) ;
        % have to interpolate in each dimension :(
        for dim = 1:3
            chordRotMat(nan_idx,dim) = interp1(x(~nan_idx), ...
                chordRotMat(~nan_idx,dim), x(nan_idx), 'makima') ;
        end
        
        % ensure normalization
        chordRotMat = chordRotMat./myNorm(chordRotMat) ;
        
        % add new entries back into full rotated and unrotated lists
        nan_ind = find(nan_idx) ;
        for q = 1:length(nan_ind)
            % get current frame, chord, and rotation matrix
            ind_curr = frame_ind(nan_ind(q)) ;
            vec_curr = chordRotMat(nan_ind(q),:) ;
            invRotMat = squeeze(body2WingArray(ind_curr,:,:))' ;
            
            % fill lists with interpolated values
            newChordRot(ind_curr,:) = vec_curr ;
            newChord(ind_curr,:) = invRotMat*vec_curr' ;
            
            % update error flag (if we turned this into nan intentionally,
            % then interpolating should solve our problem. if it was nan to
            % begin with, then the rotation matrix is likely bad, so still
            % list the error)
            error_flag(ind_curr) = ~error_flag(ind_curr) ;
        end
        
        % also interp through diagonals
        newDiag1(frame_ind(nan_idx)) = interp1(frame_ind(~nan_idx), ...
            newDiag1(frame_ind(~nan_idx)), frame_ind(nan_idx),'makima') ;
    end
    % -------------------------------------------------------------------
    % second loop -- check for irregularities in rotated frame
    % UNDER CONSTRUCTION
    %     dotDiff = [1 ; dot(newChordRot(frame_ind(1:end-1),:), ...
    %         newChordRot(frame_ind(2:end),:), 2)] ;
    %     for it = frame_ind'
    %
    %     end
    
end

% update mean and std of wing diagonals
ind = (newDiag1>0) ; % values of 0 correspond to cases where wing was not identified...
meanDiag = mean(newDiag1(ind)) ;
stdDiag  = std(newDiag1(ind)) ;
% -------------------------------------------------------------------------
%% first pass on frames with the SLOWER wing-tip.
% as above, first make some basic checks to see if we need to recalculate
% the chord and altChord. Then, in the second and third passes, we'll check
% more detailed heuristics and continuity, respectively

% find the slow intervals
slow_bouts = idx_by_thresh(~fast_ind) ;
N_slow_bouts = length(slow_bouts) ;
tflip_list = nan(N_slow_bouts,1) ;

% fit splines for tip velocity x, y, z (to be used to finely-sample)
c_vx = fit(tvec, tipVel(:,1), 'smoothingspline') ;
c_vy = fit(tvec, tipVel(:,2), 'smoothingspline') ;
c_vz = fit(tvec, tipVel(:,3), 'smoothingspline') ;

% ---------------------------------
% loop over slow bouts
for k = 1 : N_slow_bouts
    t1 = slow_bouts{k}(1) ; % interval start time (in frames)
    t2 = slow_bouts{k}(end) ; % interval end time
    v1 = tipVel_plane(t1,:) ; % tip velocities in start/end times
    v2 = tipVel_plane(t2,:) ;
    
    tvec_fine = (tvec(t1): (dt/100) : tvec(t2))' ;
    
    % find mean stroke plane normal vector - n
    if (t2>t1)
        n = mean( strokePlaneNormals(t1:t2,:)) ;
    else % t1==t2
        disp(['note t1==t2==' num2str(t1)]) ;
        error_flag(t1) = true ;
        n = strokePlaneNormals(t1,:) ;
    end
    
    n = n / norm(n) ;
    
    % evaluate tip velocity at more closely spaced time points
    tip_vx_fine = c_vx(tvec_fine) ;
    tip_vy_fine = c_vy(tvec_fine) ;
    tip_vz_fine = c_vz(tvec_fine) ;
    
    % find the wing-tip velocity component in the mean stroke-plane (fine time intervals).
    tip_v_fine       = [ tip_vx_fine, tip_vy_fine, tip_vz_fine] ;
    nmat             = repmat (n, length(tvec_fine), 1) ;
    try
        tip_v_fine_plane = tip_v_fine - ...
            nmat .* repmat ( dot(tip_v_fine, nmat, 2) , 1 , 3) ;
    catch
        keyboard ;
    end
    speed_fine_plane = myNorm(tip_v_fine_plane) ;
    
    [~, indmin] = min(speed_fine_plane) ;
    
    % the time at which the minimum speed occurs. assume this is the flip time (!)
    tflip = tvec_fine(indmin) ;
    tflip_list(k) = tflip ;
    
    % -------------------------------------------------------------------
    % go over the times in the current interval and decide whether or not
    % to recalculate chord based on velocity and diagonal length
    for it = t1:t2
        if (tvec(it)<=tflip)
            v_use = v1 ;
        else
            v_use = v2 ;
        end
        % the chord should align more with v_use
        
        c1 = mainChord(it,:) ;
        c2 = altChord(it,:) ;
        
        dot1 = dot(c1, v_use) ;
        
        % check alignment with velocity and diagonal size
        if (dot1 < 0.0) || ((diag1(it)-meanDiag)/stdDiag < -1)
            swap_flag(it) = true ;
        end
        
        % --------------------------------------------------------------
        % if true, it means there's something off about the chords, so
        % re-calculate and try to get the right ones
        if (swap_flag(it))
            if (debugFlag1)
                h_chord = showChordAndVox(data_bodyFrame, wingSide, it, c1, ...
                    c2, span(it,:), (1/voxelSize)*wingTip(it,:), ...
                    tipVel_plane(it,:))  ;
                keyboard
            end
            
            [chord_new, altChord_new, diag1_new, diag2_new, error_flag(it)] = ...
                reCalcChord(data_bodyFrame, wingSide, it, c1, c2, ...
                span(it,:), v_use, n, 'slow') ;
            
            swap_flag(it)     = true ;
            newChord(it,:)    = chord_new ; % swap
            newAltChord(it,:) = altChord_new ;
            newDiag1(it)      = diag1_new ;
            newDiag2(it)      = diag2_new ;
        else
            newChord(it,:)    = c1 ; % do not swap
            newAltChord(it,:) = c2 ;
            newDiag1(it)      = diag1(it) ;
            newDiag2(it)      = diag2(it) ;
        end
    end
end

% update mean and std of wing diagonals
ind = (newDiag1>0) ; % values of 0 correspond to cases where wing was not identified...
meanDiag = mean(newDiag1(ind)) ;
stdDiag  = std(newDiag1(ind)) ;

% -------------------------------------------------------------------------
%% second pass on SLOW FRAMES -- check heuristics

% loop over slow bouts again
for k = 1 : N_slow_bouts
    t1 = slow_bouts{k}(1) ; % interval start time (in frames)
    t2 = slow_bouts{k}(end) ; % interval end time
    v1 = tipVel_plane(t1,:) ; % tip velocities in start/end times
    v2 = tipVel_plane(t2,:) ;
    tflip = tflip_list(k) ; % flip time
    
    % ---------------------------------------------------------------------
    % go over the times in the current interval and decide if to swap the
    % chord or not. First check diagonal ratios, then see which aligns best
    % with velocity
    for it = t1:t2
        if (tvec(it)<=tflip)
            v_use = v1 ;
        else
            v_use = v2 ;
        end
        % the chord should align more with v_use
        
        c1 = mainChord(it,:) ;
        c2 = altChord(it,:) ;
        dot1 = dot(c1, v_use) ;
        dot2 = dot(c2, v_use) ;
        
        continue_process = true ;
        if ( abs(tvec(it)-tflip)*data.params.fps <= 0.25) % if very close to the flipping point <0.25 frames
            % take the longer diag
            if (diag2(it)/diag1(it)>1.2)
                swap_flag(it) = true ;
            end
            continue_process = false ;
        end
        
        if (continue_process)
            % if diag criterion applies, swap anyway, ignoring velocity criterion
            if (diag1(it) < 0.66*diag2(it)) &&  ... % if diag1 is smaller than half diag2
                    abs(diag2(it)-meanDiag)/stdDiag < 3   % and if diag2 is with 3 sigma from mean diag
                swap_flag(it) = true ;
                continue_process = false ;
            end
        end
        
        if (continue_process)
            % if diag criterion applies, swap anyway, ignoring velocity criterion
            if (diag2(it) < 0.66*diag1(it)) &&  ... % if diag1 is smaller than half diag2
                    abs(diag1(it)-meanDiag)/stdDiag < 3   % and if diag2 is with 3 sigma from mean diag
                swap_flag(it) = false ;
                continue_process = false ;
            end
        end
        
        if (continue_process)
            if (dot2>0 && dot1<0)
                swap_flag(it) = true;
            end
        end
        
        % -------------------------------------
        % swap or not
        if (swap_flag(it))
            swap_flag(it)     = true ;
            newChord(it,:)    = c2 ; % swap
            newAltChord(it,:) = c1 ;
            newDiag1(it)      = diag2(it) ;
            newDiag2(it)      = diag1(it) ;
        else
            newChord(it,:)    = c1 ; % do not swap
            newAltChord(it,:) = c2 ;
            newDiag1(it)      = diag1(it) ;
            newDiag2(it)      = diag2(it) ;
        end
    end
end

% -------------------------------------------------------------------------
%% third pass on SLOW FRAMES -- check continuity
% loop over slow bouts again. check dot product between subsequent frames.
% should be small. the cross product between subsequent chords should have
% consistent direction, but should flip at the point of wing flipping
% figure ; 
% hold on
% for k = 1 : N_slow_bouts
%     chordMat = newChord(slow_bouts{k},:) ; % chord vectors in current bout
%     dotDiff = [1 ; dot(chordMat(1:end-1,:), chordMat(2:end,:),2) ] ;
%     x = linspace(0, 1, length(dotDiff)) ; 
%     plot(x, dotDiff)
%     disp(k)
% end

% -------------------------------------------------------------------------
%% final pass on ALL FRAMES -- check "pitch angle"
% calculate something like pitch angle (the angle the chord makes with the
% stroke plane. not acutally how we do it, but close enough)
dot1 = dot(newChord, strokePlaneNormals, 2) ;
ang1 = acos(dot1) * 180 / pi;

% find any frames that have nan entries for chord, threw errors, or are
% picked up as outliers by a hampel filter
nan_idx = isnan(ang1) ; 
[~, hampel_idx] = hampel(ang1) ; 

interp_idx = nan_idx | hampel_idx | error_flag ; 

% interpolate vectors and diagonals through bad frames
x = (1:N)' ; 
for dim = 1:3
    newChord(interp_idx,dim) = interp1(x(~interp_idx), ...
        newChord(~interp_idx,dim), x(interp_idx), 'makima') ;
    newAltChord(interp_idx,dim) = interp1(x(~interp_idx), ...
        newAltChord(~interp_idx,dim), x(interp_idx), 'makima') ;
end

newDiag1(interp_idx) = interp1(x(~interp_idx), newDiag1(~interp_idx), ...
    x(interp_idx), 'makima') ;
newDiag2(interp_idx) = interp1(x(~interp_idx), newDiag2(~interp_idx), ...
    x(interp_idx), 'makima') ;

% -------------------------------------------------------------------
%% ensure that new chord is perpendicular to span
newChord = newChord - dot(newChord, span,2).*span ;
newChord = newChord./myNorm(newChord) ;
newAltChord = newAltChord - dot(newAltChord, span,2).*span ;
newAltChord = newAltChord./myNorm(newAltChord) ;

% -------------------------------------------------------------------
%% rotate back to lab frame 
% to do this, add to data structure and then transform data back to lab
% frame
data_bodyFrame.([wingSide 'ChordHats'])  = newChord ;
data_bodyFrame.([wingSide 'ChordAltHats'])  = newAltChord ;
data_bodyFrame.(['chord' wingSideNum 'AltHats']) = newAltChord ;
data_bodyFrame.(['diag' wingSideNum '1' upper(wingSide(1)) ...
    wingSide(2:end)]) = newDiag1 ;
data_bodyFrame.(['diag' wingSideNum '2' upper(wingSide(1)) ...
    wingSide(2:end)]) = newDiag2 ;

data_out = bodyToLabFrame(data_bodyFrame) ;
% -------------------------------------------------------------------
%% PLOTS (need to update!)
if (plotFlag)
    figure;
    subplot(2,3,1) ; hold on ;
    xx = 1:N ; % tvec_ms ;
    
    plot(xx, wingTip(:,1),'bo') ;
%     plot(xx, tipx_smooth,'k-','linewidth',2) ;
%     plot(xx, tipx_smooth+ESTERR,'k--') ;
%     plot(xx, tipx_smooth-ESTERR,'k--') ;
    title('wing tip x') ;
    
    subplot(2,3,2) ; hold on ;
    plot(xx, wingTip(:,2),'bo') ;
%     plot(xx, tipy_smooth,'k-','linewidth',2) ;
%     plot(xx, tipy_smooth+ESTERR,'k--') ;
%     plot(xx, tipy_smooth-ESTERR,'k--') ;
    title('wing tip y') ;
    
    subplot(2,3,3) ; hold on ;
    plot(xx, wingTip(:,3),'bo') ;
%     plot(xx, tipz_smooth,'k-','linewidth',2) ;
%     plot(xx, tipz_smooth+ESTERR,'k--') ;
%     plot(xx, tipz_smooth-ESTERR,'k--') ;
    title('wing tip z') ;
    
    subplot(2,3,4) ;
    plot(xx, tipVel(:,1),'r-','linewidth',2) ;
    title('wing tip v_x') ;
    
    subplot(2,3,5) ;
    plot(xx, tipVel(:,2),'r-','linewidth',2) ;
    title('wing tip v_y') ;
    
    subplot(2,3,6) ;
    plot(xx, tipVel(:,3),'r-','linewidth',2) ;
    title('wing tip v_z') ;
    
    
    figure ;
    plot(xx, speed,'m-^','linewidth',2) ;
    title('wing-tip speed w.r.t. body CM') ;
    ylabel('speed [m/s]') ; xlabel('time [ms]') ;
    grid on ;
    
    % generate a plot for every wing-stroke (groups of 2*DELTA frames)
    startInd = 1:(2*DELTA):N ;
    endInd   = [startInd(2:end), N] ;
    for f=1:length(startInd) 
        t1 = startInd(f) ;
        t2 = endInd(f) ;
        ind = t1:t2 ;
        figure; hold on ;
        plot3(wingTip(ind,1), wingTip(ind,2), wingTip(ind,3),'bo-') ;
        
        plot3(wingTipOther(ind,1), wingTipOther(ind,2), ...
            wingTipOther(ind,3) ,'ko') ;
        
        %plot3(tipx_smooth(ind), tipy_smooth(ind), tipz_smooth(ind),'b.-') ;
        plot3(wingTip(t1,1), wingTip(t1,2), wingTip(t1,3),'ks',...
            'markerfacecolor','r') ;
        
        a = data.AHat(t1,:) * 40 * data.params.voxelSize ;
        
        plot3( [1 -1]*a(1) , [1 -1]*a(2) , [1 -1]*a(3) ,'-','linewidth',3,'color',[0 0.7 0]) ;
        plot3(0,0,0,'o','color',[0 0.7 0],'linewidth',3) ;
        C = 0.7e-3 ;
        for it=t1:t2
            col = [0, 0, 0]*swap_flag(it) + [1 0 0]*(~swap_flag(it)) ;
            if (~fast_ind(it) && ~swap_flag(it))
                col = [1 0.5 0] ;
            end
            
            plot3( wingTip(it,1)+[0, newChord(it,1)]*C, ...
                wingTip(it,2)+[0, newChord(it,2)]*C, ...
                wingTip(it,3)+[0, newChord(it,3)]*C, '-','linewidth',2,'color',col) ;
            text( wingTip(it,1), wingTip(it,2), wingTip(it,3),['  ' num2str(it)]) ;
        end
        
        axis equal ; box on ; grid on ; view(-51, 22) ; rotate3d on ;
        axis tight;
        title(['Frame indices ' num2str(t1) ' - ' num2str(t2)]) ;
    end
    
end


%keyboard ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%% Fit affine plane to points
% this function was taken from the MathWorks website
% http://www.mathworks.com/matlabcentral/fileexchange/43305-plane-fit/content/affine_fit.m
function [n,V,p] = affine_fit(X)
%Computes the plane that fits best (lest square of the normal distance
%to the plane) a set of sample points.
%INPUTS:
%
%X: a N by 3 matrix where each line is a sample point
%
%OUTPUTS:
%
%n : a unit (column) vector normal to the plane
%V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
%plane
%p : a point belonging to the plane
%
%NB: this code actually works in any dimension (2,3,4,...)
%Author: Adrien Leygue
%Date: August 30 2013

%the mean of the samples belongs to the plane
p = mean(X,1);

%The samples are reduced:
R = bsxfun(@minus,X,p);
%Computation of the principal directions if the samples cloud
[V,D] = eig(R'*R);
%Extract the output from the eigenvectors
n = V(:,1);
V = V(:,2:end);
end

% -------------------------------------------------------------------------
%% visualize chord and alt chord
function h_chord = showChordAndVox(data, wingSide, frameNum, chord, ...
    chordAlt, span, tip, tipVel)
% ------------------------------------------------------------------
% get wing voxels
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;

row_start = frameStartInd(frameNum) ;
row_end = frameEndInd(frameNum) ;

% load in all voxel coordinates and their body part indices
coords = data.res(row_start:row_end,2:4) ;
IDX = data.RESIDX(row_start:row_end,:) ;
% get coord rows corresponding to left and right wings
wingRows = (IDX(:,data.([wingSide 'WingInd']))==1) ;

% get wing coords (in voxels, so int16)
wingVox = coords(wingRows, :) ;
cm = data.([wingSide 'WingCM'])(frameNum,:) ;

% normalize tip velocity vector
tipVelNorm = tipVel./myNorm(tipVel) ;
% voxel size
voxelSize = data.params.voxelSize ;
% --------------------------------------------------------------------
% make plot
scale = 4 ; % 2; % 4

if strcmp(wingSide,'right')
    voxColor = 'r' ;
    chordColor = 'b' ;
    chordAltColor = 'm' ;
elseif strcmp(wingSide,'left')
    voxColor = 'b' ;
    chordColor = 'r' ;
    chordAltColor = 'c' ;
else
    keyboard
end

% initialize figure window
h_chord = figure ;
hold on

% plot voxels
plot3(wingVox(:,1), wingVox(:,2), wingVox(:,3), '.','MarkerSize',3,...
    'Color',voxColor) ;

% plot span
line([cm(1), cm(1)+ scale*10*span(1)], ...
    [cm(2), cm(2)+ scale*10*span(2)], ...
    [cm(3), cm(3)+ scale*10*span(3)],...
    'Color','k','LineWidth',4);

% wing vectors (chord)
line([cm(1), cm(1)+ scale*8*chord(1)],...
    [cm(2), cm(2)+ scale*8*chord(2)],...
    [cm(3), cm(3)+ scale*8*chord(3)],...
    'Color',chordColor,'LineWidth',4);

% wing vectors (alt chord)
line([cm(1), cm(1)+ scale*8*chordAlt(1)],...
    [cm(2), cm(2)+ scale*8*chordAlt(2)],...
    [cm(3), cm(3)+ scale*8*chordAlt(3)],...
    'Color',chordAltColor,'LineWidth',4,'LineStyle','--');

% plot tip velocity
line(cm(1) + [0, scale*8*tipVelNorm(1)],...
    cm(2) + [0, scale*8*tipVelNorm(2)],...
    cm(3) + [0, scale*8*tipVelNorm(3)],...
    'Color',0.5*[1 1 1],'LineWidth',2,'LineStyle','-');

% axis properties
ax = gca ;
box(ax, 'on')
grid(ax, 'on')
axis(ax, 'equal')
rotate3d on


end

% -------------------------------------------------------------------------
%% re-calculate chord
function [chord_out, chordAlt_out, diag1, diag2, errorFlag] = ...
    reCalcChord(data, wingSide, frameNum, chord,  chordAlt, span, ...
     tipVel, strokePlaneNormal, frameType)
% ---------------------------------------
%% inputs
if ~exist('strokePlaneNormal', 'var') || isempty(strokePlaneNormal)
    strokePlaneNormal = [0, 0, 1] ;
end
if ~exist('frameType', 'var') || isempty(frameType)
    frameType = 'fast' ; % 'fast' | 'slow' (determines how we pick swap criteria)
end

% debug flags:
debugFlag1 = false ; % plot 3D voxels
debugFlag2 = false ;  % plot wing cross section with diagonals

% ---------------------------------------
%% misc calculations and params
tipVelNorm =  tipVel./myNorm(tipVel) ; % normalize wing tip velocity
delta  = 2.0;  % strip width used in finding the wing chord
N_b_pts = 200 ; % number of points to use in resampling boundary
%boundary_s = 0.5 ; % 0.1 ; % convexity of boundary calculation (0=convex hull, 1.0=super tight to points)
fitType = 'smoothingspline' ; % fit type for boundary
inContourThresh = 0.9 ; % frac. points required inside contour boundary
diagDotThresh = 0.5 ; % threshold of dot prod. for "different" wing diags
minObjSize = 100 ; %-1; % for cleaning voxels, how large must a CC be to keep

% initialize error flag
errorFlag = false ; 

% ------------------------------------------------------------------
%% get wing voxels
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;

row_start = frameStartInd(frameNum) ;
row_end = frameEndInd(frameNum) ;

% load in all voxel coordinates and their body part indices
coords = data.res(row_start:row_end,2:4) ;
IDX = data.RESIDX(row_start:row_end,:) ;
% get coord rows corresponding to left and right wings
wingRows = (IDX(:,data.([wingSide 'WingInd']))==1) ;

% get wing coords (in voxels, so int16)
wingVox = coords(wingRows, :) ;
Nvox = size(wingVox,1) ;
Rcm = data.([wingSide 'WingCM'])(frameNum,:) ;

clear df

% ----------------------------------------------------------------------
%% check that we have sufficient total wing voxels and a value for span 
if (Nvox < 100) || any(isnan(span))
    fprintf('Error with wing voxels and/or span \n')
    chord_out = nan(1,3) ;
    chordAlt_out = nan(1,3) ;
    diag1 = nan ;
    diag2 = nan ;
    errorFlag = true ;
    return
end
% ------------------------------------------------------
%% transform wing voxels
%[~, ~, keep_idx] = farthestPoint(wingVox, [0, 0, 0], LL) ;
%wingVoxCleaned = removeHullGlobs(wingVox(keep_idx,:), 400) ;
wingVoxCleaned = wingVox ;

% center the voxels
%Rcm = mean(wingVoxCleaned) ;
wingVoxTrans = wingVoxCleaned - int16(Rcm) ;

% take a strip of voxels in the wing cross section that are within some
% distance of the origin, along the span direction
spanwiseDist = dot(double(wingVoxTrans), repmat(span, Nvox,1),2) ;
mid_strip_idx = (abs(spanwiseDist) < delta) ;

% make sure we're getting points back
if (sum(mid_strip_idx) < 10)
    % first try a larger delta
    mid_strip_idx = (abs(wingVoxTrans(:,1)) < 3*delta) ;
    if (sum(mid_strip_idx) < 10)
        % if it's still a no-go, through error message and return
        fprintf('Insufficient number of voxels in strip \n')
        chord_out = nan(1,3) ;
        chordAlt_out = nan(1,3) ;
        diag1 = nan ;
        diag2 = nan ;
        errorFlag = true ; 
        return
    end
end

% just take voxels from strip, then remove any extra bits
stripVox = wingVoxTrans(mid_strip_idx,:) ;
stripVox = removeHullGlobs(stripVox, minObjSize) ;

% -----------------------------------------------------------
%% rotate voxels into (unpitched) wing frame
phi = atan2(span(2), span(1)) ;
theta = asin(span(3)) ;
rotM = body2WingRotMat(phi, theta, 0) ;

stripVox = (rotM*double(stripVox)')' ;

% also rotate other vectors
tipVelNormRot = (rotM*tipVelNorm')' ;
chordRot = (rotM*chord')' ;
chordAltRot = (rotM*chordAlt')' ;
spanRot = (rotM*span')' ;
strokePlaneNormalRot = (rotM*strokePlaneNormal')' ;

% -------------------------------------------------
% visualize in 3D?
if (debugFlag1)
    figure ;
    scale = 4 ;
    hold on
    % plot voxels
    wingVoxTransRot = (rotM*double(wingVoxTrans)')' ;
    plot3(wingVoxTransRot(:,1), wingVoxTransRot(:,2), wingVoxTransRot(:,3),...
        '.', 'MarkerSize',3, 'Color',0.6*[1 1 1]) ;
    plot3(stripVox(:,1), stripVox(:,2), stripVox(:,3), '.', ...
        'MarkerSize',5, 'Color','g') ;
    % plot span
    line([0, scale*8*spanRot(1)],...
        [0, scale*8*spanRot(2)],...
        [0, scale*8*spanRot(3)],...
        'Color','b','LineWidth',2,'LineStyle','-');
    
    % plot chord
    line([0, scale*8*chordRot(1)],...
        [0, scale*8*chordRot(2)],...
        [0, scale*8*chordRot(3)],...
        'Color','k','LineWidth',2,'LineStyle','-');
    
    % plot alt chord
    line([0, scale*8*chordAltRot(1)],...
        [0, scale*8*chordAltRot(2)],...
        [0, scale*8*chordAltRot(3)],...
        'Color','k','LineWidth',2,'LineStyle','--');
    
    % plot tip vel
    line([0, scale*8*tipVelNormRot(1)],...
        [0, scale*8*tipVelNormRot(2)],...
        [0, scale*8*tipVelNormRot(3)],...
        'Color','r','LineWidth',2,'LineStyle','-.');
    
    % axis properties
    ax = gca ;
    box(ax, 'on'); grid(ax, 'on') ; axis(ax, 'equal')
    rotate3d on
    xlabel('X') ; ylabel('Y') ; zlabel('Z')
    keyboard
end

% ----------------------------------------------
%% fit 2D boundary to projected points
% project strip voxels into 2D
stripVoxProj = [stripVox(:,2), stripVox(:,3)] ;
Nvox = size(stripVoxProj,1) ;
% use ksdensity on projected voxel points to 1) remove effect of noise and
% 2) extract boundary
[fks, xi] = ksdensity(stripVoxProj) ;
N_ks_pts = length(fks) ;
ks_mat_size = round(sqrt(N_ks_pts)) ;
fks_mat = reshape(fks, ks_mat_size, ks_mat_size) ;
y_ks = unique(xi(:,1)) ;
z_ks = unique(xi(:,2)) ;

% get contour matrix for this kernel density estimate--we'll use one of the
% countours eventually as the boundary of the projected voxels
contour_mat = contourc(y_ks, z_ks, fks_mat) ;
[vertex_list, level_list, ~] = myGetContourLines(contour_mat) ;
% check contour line levels for repeated values -- this corresponds to
% potentially disconnected contours at the same probability level. we want
% to take the one that contains the most points (~ largest)
[unique_levels, ~, ic] = unique(level_list) ;
N_unique = length(unique_levels) ; 
if (N_unique < length(level_list))
    vertex_list_sizes = cellfun(@(y) size(y,1), vertex_list) ;
    vertex_list_new = cell(1, N_unique) ;
    for q = 1:N_unique
        lvl_ind = find(q == ic) ;
        [~, max_ind] = max(vertex_list_sizes(lvl_ind)) ;
        vertex_list_new{q} = vertex_list{lvl_ind(max_ind)} ;
    end
    vertex_list = vertex_list_new ;
    %level_list = unique_levels ;
end

% get contour that contains a sufficient number of points
in_frac = 0 ;
cont_cc = length(vertex_list) + 1 ;
while (in_frac < inContourThresh) && (cont_cc > 1)
    % get vertices for next largest contour
    cont_cc = cont_cc - 1 ;
    try
        verty = vertex_list{cont_cc}(:,1) ;   
    catch
        keyboard
    end
    vertz = vertex_list{cont_cc}(:,2) ;
    
    % see how many projected voxels are contained within current contour
    in_flag = inpolygon(stripVoxProj(:,1), stripVoxProj(:,2), verty, vertz) ;
    in_frac = sum(in_flag)/Nvox ;
end

% verty and vertz should now be the boundary for our points. fit to a
% spline so we can sample more finely
curve_param = linspace(0,1, length(verty)+1)' ;
c_y = fit(curve_param, [verty ; verty(1)], fitType) ;
c_z = fit(curve_param, [vertz ; vertz(1)], fitType) ;

curve_param_fine = linspace(0,1,N_b_pts) ;
y_bound_smooth = c_y(curve_param_fine) ;
z_bound_smooth = c_z(curve_param_fine) ;

boundary_pts = [y_bound_smooth, z_bound_smooth] ;

% ----------------------------------------------------------------------
%% create a list of pairs of points ordered by their pairwise distance
% (order by largest to smallest distance)
pairDistMat = squareform(pdist(boundary_pts,'mahalanobis')) ;
%pairDistMat = squareform(pdist(boundary_pts)) ;
N_pairs = floor(size(pairDistMat,1)/2) ;
pair_vec_list = nan(N_pairs,2) ;
pair_ind_list = nan(N_pairs,2) ;
dist_list = nan(N_pairs,1) ;
for k = 1:N_pairs
    % find current maximum pairwise distance
    [~, max_ind] = nanmax(pairDistMat,[],'all','linear') ;
    [row_curr, col_curr] = ind2sub(size(pairDistMat),max_ind) ;
    % sort to make things consistent across pairs
    pair_ind = sort([row_curr, col_curr]) ;
    
    % store points and distance
    pair_ind_list(k,:) = pair_ind ;
    pair_vec_list(k,:) = boundary_pts(pair_ind(1),:) - ...
        boundary_pts(pair_ind(2),:) ;
    dist_list(k) = norm(pair_vec_list(k,:)) ;
    
    % block off those points for next round
    pairDistMat(row_curr,:) = nan ;
    pairDistMat(:, col_curr) = nan ;
end

% ---------------------------------------------------------------------
%% get top two "unique" pairwise difference vectors
% we're going to take the first large diagonal, then the next one that is
% substantially different. most distant pair becomes chord
chord = [0, pair_vec_list(1,:)] ;
diag1 = norm(chord) ;
chordHat = chord./diag1 ;

% make sure chord is pointed in positive z direction
if (dot(chordHat, strokePlaneNormalRot) < 0.0)
    chordHat = -1*chordHat ;
end
% make sure chord is pointed in direction of velocity
% if (dot(chordHat, tipVelNormRot) < 0.0)
%     chordHat(2) = -1*chordHat(2) ;
% end
% -----------------------------------------------------------------------
% check for other long diagonals that point in different direction from
% chord
dotDiff = dot(repmat(pair_vec_list(1,:),N_pairs-1,1), ...
    pair_vec_list(2:end,:),2) ;
dotDiff = dotDiff./(dist_list(1).*dist_list(2:end)) ; % normalize by vector magnitude

switch_idx = find(abs(dotDiff) < diagDotThresh,1,'first') + 1 ;

% next distant pair that is different becomes alt chord. if we find a long
% diagonal pointing in different direction, use that. otherwise take the
% diagonal that is most different
if isempty(switch_idx)
    [~, switch_idx] = min(dotDiff) ;
end
chordAlt = [0, pair_vec_list(switch_idx,:)] ;
diag2 = norm(chordAlt) ;
chordAltHat = chordAlt./diag2 ;

% make sure alt chord is pointed in positive z direction
if (dot(chordAltHat, strokePlaneNormalRot) < 0.0)
    chordAltHat = -1*chordAltHat ;
end

% ---------------------------------------------------------------------
%% check if we need to swap chord and alt chord
% we'll do this based on tip velocity in "fast" frames and more on diagonal
% in "slow" frames
velCheck1 = dot(chordHat, tipVelNormRot) ;
velCheck2 = dot(chordAltHat, tipVelNormRot) ;
swapFlag = false ;
switch frameType
    case 'fast'
        % if the alt chord lines up with velocity, but the main chord
        % doesn't, swap them. (also make sure diagonals make sense,
        % otherwise we end up with  chords pointing almost exactly in
        % direction of velocity)
        if (velCheck1 < velCheck2) && (diag2 > 1.25*diag1) || ...
                ((velCheck1 < 0.0) && (velCheck2 > 0.0))
            swapFlag = true ;
        end
    case 'slow'
        if ((diag1 < diag2) && (velCheck2 > 0)) || ...
                ((velCheck1 < 0.0) && (velCheck2 > 0.0))
            swapFlag = true ;
        end
    otherwise
        fprintf('Invalid frame type: %s \n', frameType)
        keyboard
end
% PERFORM SWAP:
if swapFlag
    % swap vector
    tmp = chordHat ;
    chordHat = chordAltHat ;
    chordAltHat = tmp ;
    
    % swap diagonal length
    tmp2 = diag1 ;
    diag1 = diag2 ;
    diag2 = tmp2 ;
end
% ---------------------------------------------
% plot to check result?
if (debugFlag2)
    figure ;
    scale = 8 ;
    hold on
    plot(stripVoxProj(:,1), stripVoxProj(:,2),'k.','markersize',5)
    %    plot(stripVoxProj(bound_ind,1), stripVoxProj(bound_ind,2), '.-',...
    %        'Color',[0, 0.6, 0])
    plot(y_bound_smooth, z_bound_smooth, 'g.-')
    % chord
    %    plot([pt2(1) pt1(1)] , [pt2(2) pt1(2)], 'r-')
    %    plot(pt1(1), pt1(2), 'ro','markerfacecolor','r')
    plot(scale*[0 chordHat(2)] , scale*[0 chordHat(3)], 'r-')
    plot(scale*chordHat(2) , scale*chordHat(3), 'ro','markerfacecolor','r')
    plot(scale*[0 chordRot(2)] , scale*[0 chordRot(3)], 'm-')
    % alt chord
    %    plot([alt_pt2(1) alt_pt1(1)] , [alt_pt2(2) alt_pt1(2)], 'b--')
    %    plot(alt_pt1(1), alt_pt1(2), 'bo','markerfacecolor','b')
    plot(scale*[0 chordAltHat(2)] , scale*[0 chordAltHat(3)], 'b--')
    plot(scale*chordAltHat(2) , scale*chordAltHat(3), 'bo','markerfacecolor','b')
    plot(scale*[0 chordAltRot(2)] , scale*[0 chordAltRot(3)], 'c--')
    plot(scale*[0 tipVelNormRot(2)] , scale*[0 tipVelNormRot(3)], 'ko-')
    % axis properties
    axis equal
    axis tight
    grid on
    box on
    keyboard
end

% ---------------------------------------------------------------------
%% rotate vectors back to body frame
chord_out = (rotM'*chordHat')' ;
chordAlt_out = (rotM'*chordAltHat')' ;

if any(isnan(chord_out))
    errorFlag = true ; 
end
%disp(frameNum) ;

end
