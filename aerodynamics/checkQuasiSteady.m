% -------------------------------------------------------------------------
% function to assess how well quasi-steady aerodynamic calculations capture
% motion of fly (single movie). going to average over wingbeats
% -------------------------------------------------------------------------
function [forcePred, torquePred, forceCalc, torqueCalc, wb, corrForce, ...
    corrTorque] = checkQuasiSteady(data, bodyFrameFlag, plotFlag)
% ---------------------------
%% inputs and params
if ~exist('bodyFrameFlag','var') || isempty(bodyFrameFlag)
    bodyFrameFlag = true ;
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true ;
end

debugFlag = false ; % show results of qs force calculations?
smoothFlag = true ; % smooth data for qs calc? good for raw data
smoothAnglesFlag = false ; % smooth body euler angles? not necessary here
addedMassFlag = false ; % calculate addedMassForce?
largePertFlag = true ; % large perturbation?

% general params for force calc
params      = defineQuasiSteadyParams() ;
g           = params.g ;
body_mass   = params.body_mass ;
% MOI_pitch   = params.MOI_pitch ;
% MOI_roll    = params.MOI_roll ;
% MOI_yaw     = params.MOI_yaw ;
beta_0      = params.beta_0 ;
span        = params.span ; 
% M2          = eulerRotationMatrix(0, -beta_0, 0) ;
% indices for body angles
defineConstantsScript ;

% low pass filter freq for body angles
body_angle_filt_lvl = 100 ; % Hz % was 50 on 1/13/20

% movie time
frames = (data.params.startTrackingTime) : (data.params.endTrackingTime) ;
dt = 1/data.params.fps ;
t = dt.*frames ;
N_frames = length(frames) ;

% % indices for angles in angleMat (should never have mixed up the order...)
% pitch_ind = 1 ;
% roll_ind = 2 ;
% yaw_ind = 3 ;

% -------------------------------------------------
%% run quasi steady force calculation
% first smooth wing angles
[~, smooth_anglesMat_R, ~, ~, ~ ] = smoothWingAngles(data, 'R') ;
[~, smooth_anglesMat_L, ~, ~, ~ ] = smoothWingAngles(data, 'L') ;

% replace raw angles with smoothed ones
data.anglesBodyFrame(:,[PHIR THETAR ETAR]) = smooth_anglesMat_R' ;
data.anglesBodyFrame(:,[PHIL THETAL ETAL]) = smooth_anglesMat_L' ;

[qs_struct, F_mat, T_mat] = calcQuasiSteady(data, debugFlag, smoothFlag, ...
    addedMassFlag) ;

% ----------------------------------------------------
%% get wingbeat cutoff times
% going to average over a window between back flip times
timesR = data.backFlipTimesR ;
timesL = data.backFlipTimesL ;

[timesR_out, timesL_out] = alignFlipTimes(timesR, timesL) ;
wbTimes = (timesR_out + timesL_out)./2 ;

% ----------------------------------------------------
%% get body angles and CM location
% center of mass position in meters
bodyCM = (data.params.voxelSize).*data.bodyCM ;

% smoothed body euler angles
[bodyPitch, bodyYaw, bodyRoll] = smoothBodyAngles(data,largePertFlag, ...
    body_angle_filt_lvl, body_angle_filt_lvl, body_angle_filt_lvl) ;

% combine in matrix and convert to radians
bodyYPR = (pi/180).*[bodyYaw, bodyPitch, bodyRoll] ; 

% -----------------------------------------------------
%% get rotation matrices
% either read from data or freshly calculate
if isfield(data,'rotM_YP') && isfield(data,'rotM_roll')
    rotM_YP = data.rotM_YP ; 
    rotM_roll = data.rotM_roll ; 
else
    [~, ~, ~, ~, ~, ~, ~, ~, ~, rotM_YP, rotM_roll, largePertFlag] = ...
        calcAnglesRaw_Sam(data, false,largePertFlag) ;
end
% need to reshape rotation matrix arrays (they need to be Nx3x3 but are
% calculated as 3x3xN)
[s1, s2, s3] = size(rotM_YP) ; 
if (s1 ~= length(t)) || (s1 < s2)
    rotM_YP = permute(rotM_YP,[3,1,2]) ; 
    rotM_roll = permute(rotM_roll,[3,1,2]) ;
end
% rotM_YP = [] ; 
% rotM_roll = [] ; 
% ------------------------------------------------------------------------
%% get wingbeat averaged predicted forces and torques based on body motion
% (in body frame)
if bodyFrameFlag
%     [forcePred, torquePred] = predictForceAndTorque(bodyCM, bodyYPR,...
%         wbTimes, t, params, largePertFlag, smoothAnglesFlag, debugFlag) ; 
    [forcePred, torquePred, wb] = predictForceAndTorque(bodyCM, bodyYPR, ...
        wbTimes, t, rotM_YP, rotM_roll, params, largePertFlag, ...
        smoothAnglesFlag, debugFlag) ; 
else
    fprintf('Under construction--calcQuasiSteady assumes body frame \n')
    keyboard
end

% ---------------------------------------------------------------
%% average quasi-steady model forces/torque over wingbeat 
% (predicted values should already be averaged)

N_wb = length(wbTimes) - 1 ;
wb0_ind = find(wbTimes < 0, 1, 'last') ;
wb = (1:N_wb) - wb0_ind ;

% initialize output arrays
F_wba            = zeros(N_wb,3) ;
T_wba            = zeros(N_wb,3) ;
F_wbaR           = zeros(N_wb,3) ;
F_wbaL           = zeros(N_wb,3) ;
% loop over wingbeats
for j = 1:N_wb
    % get indices for cutOffTimes in data matrices
    [~, ind1] = min(abs(t - wbTimes(j))) ;
    [~, ind2] = min(abs(t - wbTimes(j+1))) ;
    
    ind = ind1:ind2 ;
    
    % average 
    F_wba(j,:) = nanmedian(F_mat(ind,:)) ;
    T_wba(j,:) = nanmedian(T_mat(ind,:)) ;
    
    F_wbaR(j,:) = nanmean(qs_struct(1).F_tot(ind,:)) ; 
    F_wbaL(j,:) = nanmean(qs_struct(2).F_tot(ind,:)) ; 
%     F_wba(j,:) = nanmedian(F_mat(ind,:)) ;
%     T_wba(j,:) = nanmedian(T_mat(ind,:)) ;
end

% -------------------------------------------------------------------------
%% forces calculated using quasi steady model
scaleFactor1 = 3.0 ; %2.4;
scaleFactor2 = 10.0 ; %10 ;
% forceCalc = F_wba./scaleFactor1 ;
forceCalc = (F_wbaR + F_wbaL)./scaleFactor1 ;
torqueCalc = T_wba./scaleFactor2 ;
forceCalcR = F_wbaR./scaleFactor1 ;
forceCalcL = F_wbaL./scaleFactor1 ;

% smooth wing-beat averaged force data
for dim = 1:3
%     forceCalc(:,dim) = smooth(forceCalc(:,dim),0.2,'rloess') ;
%     torqueCalc(:,dim) = smooth(torqueCalc(:,dim),0.2,'rloess') ;
     forceCalc(:,dim) = smooth(forceCalc(:,dim),3) ;
    torqueCalc(:,dim) = smooth(torqueCalc(:,dim),3) ;
    forceCalcR(:,dim) = smooth(forceCalcR(:,dim),3) ;
    forceCalcL(:,dim) = smooth(forceCalcL(:,dim),3) ;
end
% -------------------------------------------------------------------------
%% plot results
if plotFlag
    %force
    figure ;
    for dim = 1:3
        subplot(3,1,dim)
        hold on
        %yyaxis left
        plot(wb, forcePred(:,dim)./(body_mass*g))
        
        %plot zero line
        plot([wb(1), wb(end)], [0, 0], 'k--','HandleVisibility','off') 
        
        % axis labels
        ylabel('F/mg')
        if dim == 3
           xlabel('Wingbeat') 
        end
        
        %yyaxis right
        plot(wb, forceCalc(:,dim)./(body_mass*g))
%         plot(wb, forceCalcR(:,dim)./(body_mass*g),'--')
%         plot(wb, forceCalcL(:,dim)./(body_mass*g),':')

        set(gca,'xlim',[-4, 10])
        %axis tight
    end
    legend({'predicted', 'quasi-steady'})
    
    % torque
    figure ;
    for dim = 1:3
        subplot(3,1,dim)
        hold on
        %yyaxis left
        plot(wb, torquePred(:,dim)./(body_mass*g*span))
        
        % plot zero line
        plot([wb(1), wb(end)], [0, 0], 'k--','HandleVisibility','off') 
        
        % axis labels
        ylabel('T/mgl')
        if dim == 3
           xlabel('Wingbeat') 
        end
        
        %yyaxis right
        plot(wb, torqueCalc(:,dim)./(body_mass*g*span))
        set(gca,'xlim',[-4, 10])
        %axis tight
        
    end
    legend({'predicted', 'quasi-steady'})
end

% ------------------------------------------------------------------------
%% get correlations between predicted and calculated forces/torques
fprintf('Force correlation: \n')
corrForce = corr(forcePred, forceCalc) 
fprintf('Torque correlation: \n')
corrTorque = corr(torquePred, torqueCalc)  
end