% -------------------------------------------------------------------------
% function to get angular velocity and acceleration from measured body
% euler angles
%
% angleMat = Nx3 matrix of body Euler angles in the form [pitch, roll, yaw]
%           %% needs to be in radians
% -------------------------------------------------------------------------
function [angleVel, angleAccel] = ...
    diffBodyEulerAngles(angleMat, dt, diffType, debugFlag, movslope_len)
% -------------------------
%% inputs
if ~exist('diffType','var') || isempty(diffType) 
    diffType = 'movingslope' ; % 'fit' | 'sgolay' | 'movingslope'
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = false ; 
end
if ~exist('movslope_len','var') || isempty(movslope_len)
    movslope_len = 71 ; %20 ; %100 ; %200
end
% -------------------------
%% params
fitType = 'poly4' ; %smoothingspline % function to fit (and then take derivatives of)

% savitzky golay filter parameters
sgolay_order = 5 ; 
sgolay_framelen = 101 ; 

% moving slope parameters
movslope_order = 2 ; 
% movslope_len = 71 ; %20 ; %100 ; %200

% ----------------------------------------------------
%% make sure angleMat has proper units and dimensions
if size(angleMat,1) < size(angleMat,2) 
    angleMat = angleMat' ; 
end

% this is just a rough check on degrees vs radians--should probably find a
% better way
meanAbs = nanmean(abs(angleMat(:,1))) ; 
if (meanAbs >= 10) 
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

% -------------------------
%% take derivative
switch diffType
    case 'fit'
        % create a fit object for each variable and use the differentiate
        % method
        for dim = 1:size(angleMat,2) 
            c_angle= fit(t, angleMat(:,dim), fitType) ;
            [angleVel(:,dim), angleAccel(:,dim)] = differentiate(c_angle, t) ;
        end
            
    case 'sgolay'
        % create savitzky-golay derivative filters
        % initialize matrix for derivatives
        derMat = zeros(size(angleMat,1),size(angleMat,2),2) ; 
        
        % get sgolay filter coeffs
        [~, g] = sgolay(sgolay_order, sgolay_framelen) ;
        for dim = 1:size(angleMat,2)
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
        for dim = 1:size(angleMat,2)
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
%     ylabels_vel = {'Pitch Vel. (rad/s)', 'Roll Vel. (rad/s)', ...
%         'Yaw Vel. (rad/s)'} ;
    for j = 1:3
        subplot(3,1,j)
        hold on
        plot(t, (1/dt)*[nan ; diff(angleMat(:,j))], '.', ...
            'Color', 0.5*[1 1 1 1])
        plot(t, angleVel(:, j), 'r-')
        
        %ylabel(ylabels_vel{j}) 
        xlabel('Time (s)')
        axis tight
    end
    
    % acceleration
    figure ;
%     ylabels_accel = {'Pitch Accel (rad/s^2)', 'Roll Accel (rad/s^2)', ...
%         'Yaw Accel (rad/s^2)'} ;
    for j = 1:3
        subplot(3,1,j)
        hold on
        plot(t, (1/dt)^2*[nan ; nan ; diff(diff(angleMat(:,j)))], '.', ...
            'Color', 0.5*[1 1 1 1])
        plot(t, angleAccel(:, j), 'r-')
        
%         ylabel(ylabels_accel{j}) 
        xlabel('Time (s)')
        axis tight
    end
    
    keyboard
end
end