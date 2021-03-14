% -------------------------------------------------------------------------
% function to remove jumps in the estimation of the body longitudinal axis
% that occur when the sign gets flipped (since the body axis is first
% estimated as the 1st principal component of the body voxels, there is a
% degeneracy in sign)
% -------------------------------------------------------------------------
function data = correctAHat(data)
% --------------------------------------
% read in body axis from data struct
AHat = data.AHat ; 
N_frames = size(AHat,1) ; 
earlyFrameWindow = 1:50 ; 
cosDistThresh = 0.02 ; % this corresponds to a little over 10 degrees (11.48)

% initialize new values for AHat
AHat_new = AHat ; 
flip_idx = false(N_frames, 1) ; 

% first see if we can apply hampel filter to catch any weird jumps
AHat_new = hampel(AHat_new, 101) ; 
AHat_new = AHat_new./myNorm(AHat_new) ; 

% check that the first 50 entries agree generally. if they do, we'll just
% march through the array and switch the signs on ones that don't match the
% previous frame. if they don't agree, then we can just try to do the same
% thing but in reverse?
earlyFrameTest = pdist(AHat_new(earlyFrameWindow,:),'cosine') ; 

if (nanmax(earlyFrameTest) < cosDistThresh)
    % in this case, the feed-forward approach should be fine
    for ind = 2:N_frames
       dotCheck = dot(AHat_new(ind-1,:), AHat_new(ind,:)) ; 
       if dotCheck < 0
          AHat_new(ind,:) = -1*AHat(ind,:) ; 
          flip_idx(ind) = true ; 
       end
    end
else
    keyboard
    % try reverse loop?
    AHat_new = flipud(AHat_new) ; 
    for ind = 2:N_frames
       dotCheck = dot(AHat_new(ind-1,:), AHat_new(ind,:)) ; 
       if dotCheck < 0
          AHat_new(ind,:) = -1*AHat(ind,:) ; 
          flip_idx(ind) = true ; 
       end
    end
    AHat_new = flipud(AHat_new) ; 
end

if (1)
   figure ; 
   hold on
   plot(asin(AHat(:,3))) 
   plot(asin(AHat_new(:,3))) 
   axis tight
   
end

% --------------------------------------------------------------
% assign new body axis vectors to struct 
data.AHat = AHat_new ; 

% in cases where we flipped AHat sign, will likely need to swap L<->R
swapFrames= find(flip_idx) ; 
data = swapWingLeftRight(data,swapFrames) ; 

end