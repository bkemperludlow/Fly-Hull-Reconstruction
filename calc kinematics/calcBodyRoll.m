%--------------------------------------------------------------------------
% function to calculate the body roll of the fly. copied and pasted from
% calcAngles_quick_and_dirty_*.m
%
%   INPUTS:
%       -rhoTimes = frame numbers wherein roll is estimated by rollVector
%       -rollVectors = Nx3 matrix of normalized vectors that point from
%       R-->L wing CM (or hinge) that are used to estimate the body roll
%       -t = time in seconds for the movie
%       -rotM_YP = 3x3xN array containing rotation matrices that will first
%       unyaw then unpitch the fly, so that its long body axis is along the
%       x axis
%       -data_params = the "params" field of the normal data strucure.
%       Contains info like the frames per second (fps), first tracked
%       frame, etc.
%
%   OUTPUTS:
%       -smoothed_rho = estimated value of body roll angle at each time
%       point in t
%       -sp_rho = spline fit used to get body roll based on the
%       measurements made at rhoTimes
%       -rho_t = times (in seconds) that correspond to the frames in
%       rhoTimes
%       -rho_samp = roll estimates at rhoTimes, used to fit spline
%
%--------------------------------------------------------------------------
function [smoothed_rho, sp_rho, rho_t, rho_samp, rotM_roll] = ...
    calcBodyRoll(rhoTimes, rollVectors, t, rotM_YP, data_params, ...
      largePertFlag)
%--------------------------------------------------------------------------
%% params
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false ; 
end
binRhoTimesFlag = true;
RAD2DEG = 180 / pi ;
DEG2RAD = pi / 180 ;
Nimages = length(t) ;
smoothingParams = setSmoothingParams() ;
rollErr = smoothingParams.roll_est_err ;
%--------------------------------------------------------------------------
%% provide endpoints for roll estimation
if (rhoTimes(1)~=1)
    rollVectors(1,:) = rollVectors(rhoTimes(1),:) ;
    rhoTimes = [1 rhoTimes] ;
end
if (rhoTimes(end)~=Nimages)
    rollVectors(end,:) = rollVectors(rhoTimes(end),:) ;
    rhoTimes = [rhoTimes Nimages] ; % when using estimateRollVectos
end

% bin together roll estimates that are close in time
if binRhoTimesFlag
    rhoTimes_diff = diff(rhoTimes) ;
    chunk_ind = find(rhoTimes_diff > 2) + 1 ;
    % in case whave high sampling or something, make sure binning doesn't
    % fuck us over
    if ~isempty(chunk_ind)
        % include endpoints in chunk indices
        chunk_ind = [1, chunk_ind, (length(rhoTimes)+1)] ;
        
        rhoTimes_old = rhoTimes ;
        rollVectors_old = rollVectors ;
        
        rhoTimes = nan(1,(length(chunk_ind)-1)) ;
        rollVectors = nan(Nimages,3) ;
        
        % loop through chunks and bin values
        for qq = 1:(length(chunk_ind)-1)
            ind1 = chunk_ind(qq) ;
            ind2 = chunk_ind(qq+1) - 1 ;
            rhoTimes(qq) = round(mean(rhoTimes_old(ind1:ind2))) ;
            rollVectors(rhoTimes(qq),:) = ...
                nanmean(rollVectors_old(rhoTimes_old(ind1:ind2),:),1) ;
        end
    end
end

%--------------------------------------------------------------------------
%% estimate roll angle at rhoTimes by undoing the yaw and pitch rotations
rho_samp = zeros(length(rhoTimes),1) ;
rotM_roll = nan(3,3, Nimages) ; 
rotRollVecs = nan(size(rollVectors)) ; 
%rho_samp2 = rho_samp ;
rho_t    = rhoTimes + data_params.firstTrackableFrame - 1 ;
rho_t    = rho_t / data_params.fps  ;

for j=1:length(rhoTimes)
    it = rhoTimes(j) ;
    rotM = squeeze(rotM_YP(:,:,it)) ;  % rotation matrix to undo yaw and pitch
    yb = rollVectors(it,:) ; %  goes from R-->L wing hinges
    
    %try to correct for wing swapping, as can sometimes occur in large perts
    if largePertFlag && (j > 1)
    	dotCheck = dot(yb, rollVectors(rhoTimes(j-1),:)) ; 
        if (dotCheck < -0.5)
            yb = -1*yb ; 
            rollVectors(it,:) = yb ;
        end
    end
    
    % perform rotation
    rotyb = rotM * yb' ;
    
    % make sure (rotated) roll vector has zero component in x direction
    rotyb(1) = 0 ; 
    rotyb = rotyb./norm(rotyb) ; 
    
    % estimate angle
    rho_samp(j) = acos( rotyb(2) ) * RAD2DEG * sign(rotyb(3)) ; 
    rotRollVecs(it,:) = rotyb' ; 
    rotM_roll(:,:,it) = eulerRotationMatrix(0,0, rho_samp(j)*DEG2RAD) ; 
end

% ----------------------------------------------
% unwrap values of rho to check for large jumps
% NB: this doesn't affect roll matrices, since they don't care about shifts
% by pi
rho_samp = RAD2DEG*unwrap(DEG2RAD*rho_samp) ;  

%--------------------------------------------------------------------------
%% fit spline through roll points
tol = rollErr ; % 4 ;
[sp_rho, ~] =  spaps(rho_t, rho_samp, tol) ;
rho0 = fnval(sp_rho, t) ;

smoothed_rho = rho0 ;

%--------------------------------------------------------------------------
%% fill in roll matrices by interpolating through roll vectors
if largePertFlag
    % if it's a large perturbation, we do a frame-by-frame wind up rotation
    % to avoid gimbal lock weirdness
    frames = (1:Nimages)' ;
    %c_x = fit(frames(rhoTimes), rotRollVecs(rhoTimes,1), 'cubicinterp') ;
    c_y = fit(frames(rhoTimes), rotRollVecs(rhoTimes,2), 'cubicinterp') ;
    c_z = fit(frames(rhoTimes), rotRollVecs(rhoTimes,3), 'cubicinterp') ;
    
    rotRollVecs_interp = [zeros(Nimages,1), c_y(frames), c_z(frames)] ;
    rotRollVecs_interp = rotRollVecs_interp ./ ...
        repmat(myNorm(rotRollVecs_interp),1,3) ;
    
    rotM_roll(:,:,1) = eulerRotationMatrix(0,0,rho_samp(1)*DEG2RAD) ;
    rotCheck = nan(Nimages,3) ;
    
    for k = 2:Nimages
        rotyb_curr =  rotRollVecs_interp(k,:) ;
        rotM_roll_prev = squeeze(rotM_roll(:,:,k-1)) ;
        
        rotyb_rot = rotM_roll_prev*rotyb_curr' ;
        rho_samp_curr = real(acos( rotyb_rot(2) )) * sign(rotyb_rot(3)) ;
        rotM_roll_curr =  eulerRotationMatrix(0,0,rho_samp_curr) ;
        rotM_roll(:,:,k) = rotM_roll_curr * rotM_roll_prev ;
        rotCheck(k,:) = rotyb_rot' ;
    end
    
    if (0)
        figure ;
        hold on
        plot(rotCheck(:,1),'-')
        plot(rotCheck(:,2),':')
        plot(rotCheck(:,3),'--')
        title('Rotated Roll Vector Check')
    end
else
    % if it's not a large perturbation, just use standard method to get
    % roll matrices
    for k = 1:Nimages
        rotM_roll(:,:,k) = eulerRotationMatrix(0,0,DEG2RAD*smoothed_rho(k)) ;
    end
end

end

%% old attemp at large pert code (from 2/4/19)
%{
if largePertFlag
    % in this case, because the fly may be flipping around, we calculate
    % roll frame to frame (or rather, rhoTime to rhoTime) to avoid big
    % skips
    rho_diff = zeros(size(rho_samp)) ;
    % perform unyaw, unpitch on first rollVector. Then calculate roll
    % angle, and rotate vector by that amount
    rollVec_init = rollVectors(rhoTimes(1),:) ;
    rotM = squeeze(rotM_YP(:,:,rhoTimes(1))) ;
    rollVec_rot1 = rotM * rollVec_init' ;
    %rho_init = acos( rollVec_rot1(2) ) * sign(rollVec_rot1(3)) ;
    rho_init = atan2(rollVec_rot1(3), rollVec_rot1(2)) ;
    rho_diff(1) = rho_init ;
    rotM_roll = eulerRotationMatrix(0,0, rho_init) ;
    rollVec_rot2 = rotM_roll * rollVec_rot1 ;
    
    % initialize some data containers
    rotM_roll_array = zeros(3,3,length(rhoTimes)) ;
    rollVec_rot_array = zeros(length(rhoTimes),3) ;
    rotM_roll_array(:,:,1) = rotM_roll ;
    rollVec_rot_array(1,:) = rollVec_rot2 ;
    
    % now loop through and get rotations between subsequent roll samplings
    for j = 2:length(rhoTimes)
        it = rhoTimes(j) ;
        rotM = squeeze(rotM_YP(:,:,it)) ;  % rotation matrix to undo yaw and pitch
        rotM_roll_prev = squeeze(rotM_roll_array(:,:,j-1)) ;
        rollVec_curr = rollVectors(it,:) ;
        
        % first do the unyaw, unpitch
        rollVec_curr_rot1 = rotM * rollVec_curr' ;
        % now perform roll rotation from previous rhoTime
        rollVec_curr_rot2 = rotM_roll_prev * rollVec_curr_rot1 ;
        
        % from this, calculate current roll angle
%         rho_curr = acos( rollVec_curr_rot2(2) ) * sign(rollVec_curr_rot2(3)) ;
        rho_curr = atan2( rollVec_curr_rot2(3), rollVec_curr_rot2(2) ) ;
        rho_diff(j) = rho_curr ;
        
        % create next rotation mat
        rotM_roll_curr = eulerRotationMatrix(0,0,rho_curr) ;
        rotM_roll_array(:,:,j) = rotM_roll_curr * rotM_roll_prev ;
        rollVec_rot_array(j,:) = (rotM_roll_curr * rollVec_curr_rot2)' ;
    end
    rho_samp = rad2deg*cumsum(rho_diff) ;
else
%}