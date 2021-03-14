% -------------------------------------------------------------------------
% quick function to convert body frame x and z coordinates to lab frame
% from simulation results
% -------------------------------------------------------------------------
function [x_lab, x_dot_lab, z_lab, z_dot_lab] = ...
    convertSimCoordsToLab(x_dot, z_dot, thetaB, dt) 
% ------------------------
% number of time points
x_dot_lab = zeros(size(x_dot)) ;
z_dot_lab = zeros(size(z_dot)) ;

% loop over each time step and rotate current velocity vector
for k = 1:length(x_dot)
    % current rotation matrix
    rotM = eulerRotationMatrix(0, thetaB(k),0) ; 
    
    % rotate velocity vector
    vel_vec = [x_dot(k) ; 0 ; z_dot(k)] ;
    vel_vec_lab = rotM'*vel_vec ;
    x_dot_lab(k) = vel_vec_lab(1) ;
    z_dot_lab(k) = vel_vec_lab(3) ;
end

% get lab frame position by integrating lab frame velocity
x_lab = dt.*cumtrapz(x_dot_lab) ; 
z_lab = dt.*cumtrapz(z_dot_lab) ; 


end