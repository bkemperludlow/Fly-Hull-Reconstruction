%--------------------------------------------------------------------------
% once we've obtained voxel shifts to align the 2 side views to the top
% view for each of the calibration wand points, we can interpolate those
% data to make an initial guess at the shift we'll need to align the fly
% voxels on any given frame
%--------------------------------------------------------------------------
function interp_vals = interpolateVoxShifts(vox_shifts, centers_3D, query_pts)

% initialize array
interp_vals = nan(size(query_pts)) ; 
%interpolants = [] ; 
% loop through shift dimensions
for i = 1:size(query_pts,2)
    F = scatteredInterpolant(centers_3D, vox_shifts(:,i)) ; 
    interp_vals(:,i) = F(query_pts) ; 
    %interpolants = [interpolants ; F] ; 
end

end