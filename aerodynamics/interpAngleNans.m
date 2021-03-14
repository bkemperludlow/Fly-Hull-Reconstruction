% -------------------------------------------------------------------------
% interpolate nan values in Nx3 angle (or other) data matrix
% -------------------------------------------------------------------------
function angleMat = interpAngleNans(angleMat)
nan_idx = any(isnan(angleMat), 2) ; 

frames = (1:size(angleMat,1))' ;
for i = 1:size(angleMat,2) 
   angleMat(:,i) = interp1(frames(~nan_idx), angleMat(~nan_idx,i), ...
       frames, 'spline') ; 
end

end