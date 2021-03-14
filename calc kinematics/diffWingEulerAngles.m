% -------------------------------------------------------------------------
% function to get angular velocity and acceleration from measured wing
% euler angles
%
% angleMat = Nx3 matrix of body Euler angles in the form 
%           [stroke, deviation, pitch]
%%% needs to be in radians
% 
% -------------------------------------------------------------------------
function [angleVel, angleAccel] = ...
    diffWingEulerAngles(angleMat, dt, diffType, debugFlag)
% -------------------------
%% inputs
if ~exist('diffType','var') || isempty(diffType) 
    diffType = 'movingslope' ; % 'fit' | 'sgolay' | 'movingslope'
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = false ; 
end

% -------------------------
%% params
fitType = 'cubicinterp' ; % 'smoothingspline' ; % function to fit (and then take derivatives of)

% savitzky golay filter parameters
sgolay_order = 3 ; 
sgolay_framelen = 11 ; 

% moving slope parameters
movslope_order = 3 ; 
movslope_len = max([round(5e-4/dt), movslope_order+1]) ; 

% ----------------------------------------------------
%% make sure angleMat has proper units and dimensions
if size(angleMat,1) < size(angleMat,2) 
    angleMat = angleMat' ; 
end

% this is just a rough check on degrees vs radians--should probably find a
% better way
strokeMax = nanmedian(abs(angleMat(:,1))) ; 
if (strokeMax >= 2*pi) 
    % if the mean of the absolute value for pitch is ~10, it's probably in
    % degrees
    angleMat = (pi/180).*angleMat ;
end

% ------------------------------
%% initialize arrays
angleVel = nan(size(angleMat)) ; 
angleAccel = nan(size(angleMat)) ; 

% also get time
t = dt*(1:size(angleMat,1))' ; 

% fill in nan values
nan_idx = any(isnan(angleMat),2) ; 
frames = (1:size(angleMat,1))' ;
for i = 1:size(angleMat,2) 
   angleMat(:,i) = interp1(frames(~nan_idx), angleMat(~nan_idx,i), ...
       frames, 'spline') ; 
end
% -------------------------
%% take derivative
switch diffType
    case 'fit'
        % create a fit object for each variable and use the differentiate
        % method
        for dim = 1:3 
           c_angle = fit(t, angleMat(:,dim), fitType) ; 
           [angleVel(:,dim), angleAccel(:,dim)] = differentiate(c_angle, t) ; 
        end
        
    case 'sgolay'
        % create savitzky-golay derivative filters
        % initialize matrix for derivatives
        derMat = zeros(size(angleMat,1),3,2) ; 
        
        % get sgolay filter coeffs
        [~, g] = sgolay(sgolay_order, sgolay_framelen) ;
        for dim = 1:3
            % read in pitch, roll, or yaw
            data_curr = angleMat(:,dim) ;
            
            % loop through levels of differentiation
            for p = 1:2
                derMat(:,dim, p) = conv(data_curr, ...
                    factorial(p)/(-dt)^p * g(:,p+1),'same');
            end
        end
        
        angleVel = derMat(:,:,1) ; 
        angleAccel = derMat(:,:,2) ; 
        
    case 'movingslope'
        % use the movingslope tool from:
        % http://www.mathworks.com/matlabcentral/fileexchange/16997-movingslope movingslope
        movslope_len = min([movslope_len, size(angleMat,1)-1]) ; 
        for dim = 1:3 
           angleVel(:,dim) = (1/dt)*movingslope(angleMat(:,dim),...
               movslope_len, movslope_order) ; 
           angleAccel(:,dim) = (1/dt)*movingslope(angleVel(:,dim),...
               movslope_len, movslope_order) ; 
        end
    
    otherwise
        fprintf('Invalid selection for diffType: %s \n', diffType)
        keyboard 
end

% ---------------------------------
%% plot results to debug?
if debugFlag
    % velocity
    figure ;
    ylabels_vel = {'Stroke Vel. (rad/s)', 'Deviation Vel. (rad/s)', ...
        'Rotation Vel. (rad/s)'} ;
    for j = 1:3
        subplot(3,1,j)
        hold on
        plot(t, (1/dt)*[nan ; diff(angleMat(:,j))], '.', ...
            'Color', 0.5*[1 1 1 1])
        plot(t, angleVel(:, j), 'r-')
        
        ylabel(ylabels_vel{j}) 
        xlabel('Time (s)')
        axis tight
    end
    
    % acceleration
    figure ;
    ylabels_accel = {'Stroke Accel (rad/s^2)', 'Deviation Accel (rad/s^2)', ...
        'Rotation Accel (rad/s^2)'} ;
    for j = 1:3
        subplot(3,1,j)
        hold on
        plot(t, (1/dt)^2*[nan ; nan ; diff(diff(angleMat(:,j)))], '.', ...
            'Color', 0.5*[1 1 1 1])
        plot(t, angleAccel(:, j), 'r-')
        
        ylabel(ylabels_accel{j}) 
        xlabel('Time (s)')
        axis tight
    end
end
end