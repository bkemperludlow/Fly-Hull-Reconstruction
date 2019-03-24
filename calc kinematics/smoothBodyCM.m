function [bodyCM_smooth, bodyVel, bodyAccel] = ...
    smoothBodyCM(bodyCM,smoothType, debugFlag)

if nargin < 2
    smoothType = 'spline' ;
    debugFlag = false ;
elseif nargin < 3
    debugFlag = false ; 
end

%% define params
FPS = 8000 ; % frame rate for videos
voxelSize = 50e-6 ; % voxel to real space conversion (50 micron per voxel side)
frames = 1:size(bodyCM,1) ; 
t = (1/FPS)*frames ; % in seconds

dt = mean(diff(t)) ; 

smoothingParams = setSmoothingParams ; 
sGolayOrder = smoothingParams.body_sGolayOrder ;    % poly order for savitzky golay
smoothWindow = smoothingParams.body_smoothWindow ; % smoothing window for savitzky golay
halfPowFreq = smoothingParams.body_cm_filt_lvl ; % half power frequency for butterworth filter
filterOrder = smoothingParams.body_cm_filt_order ;  % filter order 
spanWindow = smoothingParams.body_spanWindow; % regression span for loess
expectedErr = smoothingParams.body_cm_est_err ; % in voxels (for spline)

%% cases for different smoothing methods
switch smoothType
%-----------------------------------
    case 'lowpass'
        %% low pass buttwerworth filter for Cm data
        d1 = designfilt('lowpassiir','FilterOrder',filterOrder,'SampleRate',FPS, ...
            'HalfPowerFrequency',halfPowFreq,'DesignMethod','butter'); 

        bodyCM_smooth = filtfilt(d1,bodyCM) ; 
%-----------------------------------
    case 'sgolay'
        %% smoothing savitzky golay filter
        bodyCM_smooth = zeros(size(bodyCM)) ; 
        for i = 1:3
            bodyCM_smooth(:,i) = smooth(bodyCM(:,i),smoothWindow, 'sgolay',sGolayOrder) ;
        end
%-----------------------------------
    case 'kalman'
        %% smoothing kalman filter
        % motionModel = 'ConstantVelocity' ; 
        % initialLocation = bodyCM(1,:) ; 
        % initialEstimateError = [1 1] ; 
        % motionNoise = 20*[1 1] ; 
        % measurementNoise = 5 ; 
        % kalman = configureKalmanFilter(motionModel,initialLocation,...
        %     initialEstimateError,motionNoise,measurementNoise);

        stateTransitionModel = [1, 1, 0, 0, 0, 0 ; 0, 1, 0, 0, 0, 0 ; ...
                                0, 0, 1, 1, 0, 0 ; 0, 0, 0, 1, 0, 0 ; ...
                                0, 0, 0, 0, 1, 1 ; 0, 0, 0, 0, 0, 1] ; 

        measurementModel     = [1, 0, 0, 0, 0, 0 ; 0, 0, 1, 0, 0, 0 ; ...
                                0, 0, 0, 0, 1, 0] ;
        sigma2               = 0.25 ; 
        processNoise         = sigma2 * [dt^3/3, dt^2/2, 0, 0, 0, 0 ; ...
                                           dt^2/2, dt, 0, 0, 0, 0 ; ...
                                           0, 0, dt^3/3, dt^2/2, 0, 0 ; ....
                                           0, 0, dt^2/2, dt, 0, 0 ; ...
                                           0, 0, 0, 0, dt^3/3, dt^2/2 ; ...
                                           0, 0, 0, 0, dt^2/2, dt ] ; 
        measurementNoise     = 50 ;  %4
        stateCovariance      = 1 ; 

        kalman = vision.KalmanFilter(stateTransitionModel,measurementModel,...
            'ProcessNoise',processNoise,'MeasurementNoise',measurementNoise,...
            'StateCovariance',stateCovariance);

        v_init = bodyCM(2,:) - bodyCM(1,:) ; 
        state_init = [bodyCM(1,1), v_init(1), bodyCM(1,2), v_init(2), ...
            bodyCM(1,3), v_init(3)] ; 
        kalman.State = state_init ; 

        bodyCM_smooth = nan(size(bodyCM)) ; 
        for k = 1:size(bodyCM,1) 
            trackedPosition = predict(kalman) ; 
            trackedPosition = correct(kalman,bodyCM(k,:));
            bodyCM_smooth(k,:) = trackedPosition ; 
        end
%-----------------------------------
    case 'loess'
        %% local regression using weighted LLS and a 2nd degree polynomial model
        expectedErr = 0.01 ; 
        bodyCM_smooth = zeros(size(bodyCM)) ; 
        for i = 1:3
            temp_smooth = smooth(t,bodyCM(:,i),spanWindow, 'loess') ;
            
            %sp_smooth = fit(t', temp_smooth, 'smoothingspline') ; 
            %bodyCM_smooth(:,i) = sp_smooth(t) ; 
            
            sp_smooth = mySplineSmooth(t,temp_smooth,expectedErr) ;
            bodyCM_smooth(:,i) = fnval(sp_smooth, t) ;
            
            %bodyCM_smooth(:,i) = smooth(t,bodyCM(:,i),spanWindow, 'loess') ;
        end

%-----------------------------------
    case 'spline'
        %% smoothing spline
        
        bodyCM_smooth = zeros(size(bodyCM)) ; 
        sp_x = mySplineSmooth(t,bodyCM(:,1),expectedErr) ;
        sp_y = mySplineSmooth(t,bodyCM(:,2),expectedErr) ;
        sp_z = mySplineSmooth(t,bodyCM(:,3),expectedErr) ;

        bodyCM_smooth(:,1) = fnval(sp_x, t) ; 
        bodyCM_smooth(:,2) = fnval(sp_y, t) ; 
        bodyCM_smooth(:,3) = fnval(sp_z, t) ; 
        
%-----------------------------------
    otherwise
        disp('No smoothing method selected')
        keyboard ;
end

%% compare smoothed results and raw data
if debugFlag 
   figPosition = [665   325   961   287] ; 
   h_pos = figure('PaperPositionMode','auto','Position',figPosition) ; 
   label_cell = {'X Position','Y Position','Z Position'} ; 
   for k = 1:3 
       subplot(1,3,k)
       hold on
       plot(bodyCM_smooth(:,k),'r-')
       plot(bodyCM(:,k),'k.') 

       axis tight
       xlabel('Frame #')
       ylabel(label_cell{k})
       title(label_cell{k})
   end 
   h_res = figure('PaperPositionMode','auto') ; 
   max_err = max(abs(bodyCM(:) - bodyCM_smooth(:))) ; 
   edges = (-1*max_err):0.05:(max_err) ;
   %Nbins = 50 ; 
   for k = 1:3 
       subplot(3,1,k)
       histogram((bodyCM(:,k) - bodyCM_smooth(:,k)),edges,'normalization','pdf')

       axis tight
       xlabel([label_cell{k} 'Residuals'])
       ylabel('PDF')
       title(label_cell{k})
   end 
end

%% fit cubic interpolant to take derivatives
bodyCM_smooth = voxelSize * bodyCM_smooth ; 
x_cm = bodyCM_smooth(:,1) ; 
y_cm = bodyCM_smooth(:,2) ; 
z_cm = bodyCM_smooth(:,3) ; 

c_x = fit(t', x_cm, 'cubicinterp') ; 
c_y = fit(t', y_cm, 'cubicinterp') ; 
c_z = fit(t', z_cm, 'cubicinterp') ; 

[vx, ax] = differentiate(c_x,t) ; 
[vy, ay] = differentiate(c_y,t) ; 
[vz, az] = differentiate(c_z,t) ; 

% combine data into matrices
bodyVel = [vx , vy , vz] ; 
bodyAccel = [ax, ay, az] ; 

%% plot results?
if debugFlag 
   figure; 
   label_cell_vel = {'X Vel','Y Vel','Z Vel'} ; 
   for k = 1:3 
       subplot(3,1,k)
       hold on
       plot(t,bodyVel(:,k),'r-')
       %plot(t,[nan; diff(bodyCM(:,k))],'k.') 

       axis tight
       xlabel('Time')
       ylabel(label_cell_vel{k})
       title(label_cell_vel{k})
   end
   %=============
   figure;
   label_cell_accel = {'X Accel','Y Accel','Z Accel'} ; 
   for k = 1:3 
       subplot(3,1,k)
       hold on
       plot(t,bodyAccel(:,k),'r-')
       %plot(t,[nan ; nan ; diff(diff(bodyAccel(:,k)))],'k.') 

       axis tight
       xlabel('Time')
       ylabel(label_cell_accel{k})
       title(label_cell_accel{k})
   end 
end

end