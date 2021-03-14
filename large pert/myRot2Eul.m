% -------------------------------------------------------------------------
% attempt at a function to get euler angles from rotation matrix (to be
% used with large pert calculations). based on:
% "Computing Euler angles from a rotation matrix" by Gregory Slabaugh
% -------------------------------------------------------------------------
function [phi1, phi2, theta1, theta2, psi1, psi2] = myRot2Eul(rotM)
%check whether theta = +/- pi/2
tol = 4*eps ; 
halfPiCheck = (abs(1 - abs(rotM(1,3))) < tol) ; % if true, theta ~ +/- pi/2

% if we're not at pitch singularity, can find 2 solutions for each angle
if ~halfPiCheck
    theta1 = -1*asin(rotM(1,3)) ;
    theta2 = pi - theta1 ;
    psi1 = atan2(rotM(2,3)/cos(theta1), rotM(3,3)/cos(theta1)) ;
    psi2 = atan2(rotM(2,3)/cos(theta2), rotM(3,3)/cos(theta2)) ;
    phi1 = atan2(rotM(2,1)/cos(theta1), rotM(1,1)/cos(theta1)) ;
    phi2 = atan2(rotM(2,1)/cos(theta2), rotM(1,1)/cos(theta2)) ;
else
    % if we are at pitch singularity, infinite solutions. Try to guess theta
    % by sign of sin(theta), but give nan for psi and phi (we'll interpolate
    % these later)
    if (sign(rotM(1,3)) < 0)
        theta1 = pi/2 ; 
        theta2 = -3*pi/2 ; 
    else
        theta1 = -1*pi/2 ; 
        theta2 = 3*pi/2 ; 
    end
    psi1 = nan ; 
    psi2 = nan ; 
    phi1 = nan ; 
    phi2 = nan ; 
end

end