% -------------------------------------------------------------------------
% function to take in a single movie's data structure and calculate the
% quasi-steady forces and torques
% -------------------------------------------------------------------------
function [qs_force_struct, F_mat, T_mat] = calcQuasiSteady(data, plotFlag,...
    smoothFlag, addedMassFlag)
% -----------------------
%% inputs and params
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true ;
end
if ~exist('smoothFlag','var') || isempty(smoothFlag)
    smoothFlag = true ;
end
if ~exist('addedMassFlag','var') || isempty(addedMassFlag)
    addedMassFlag = false ;
end

debugFlag = false ; % plot within qs force struct function?

DEG2RAD = (pi/180) ; % angle conversion
%RAD2DEG = (180/pi) ;

defineConstantsScript ; % indices to access body/wing angles
params = defineQuasiSteadyParams() ; % parameters for qs calculations
% -----------------------------
%% read data from struct
% time in seconds
frames = (data.params.startTrackingTime : data.params.endTrackingTime) ;
dt = (1/data.params.fps) ;
t = dt.*frames ; % seconds

% wing angles (in radians)
wingAngRight = DEG2RAD.*(data.anglesBodyFrame(:,[PHIR, THETAR, ETAR])) ;
wingAngLeft  = DEG2RAD.*(data.anglesBodyFrame(:,[PHIL, THETAL, ETAL])) ;

% body center of mass (in meters)
bodyCM = (data.params.voxelSize) .* (data.bodyCM) ;

% body euler angles (in radians)
bodyYPR = DEG2RAD .* (data.anglesLabFrame(:,[PHIB, THETAB, RHO])) ;

% -----------------------------
%% calculate forces
qs_force_structR =quasiSteadyForceAll(wingAngRight, t, 'R', ...
    params, bodyCM, bodyYPR, debugFlag, smoothFlag, addedMassFlag) ;
qs_force_structL =quasiSteadyForceAll(wingAngLeft, t, 'L', ...
    params, bodyCM, bodyYPR, debugFlag, smoothFlag, addedMassFlag) ;

% -----------------------------
%% calculate torques
qs_force_structR = quasiSteadyTorque(qs_force_structR) ;
qs_force_structL = quasiSteadyTorque(qs_force_structL) ;

% combine data structures
qs_force_struct = vertcat(qs_force_structR, qs_force_structL) ;

% ------------------------------
%% combine (and plot?) results
if plotFlag
    % plot forces and torques
    h_qs = figure ;
    gap = [0.05, 0.08] ;
    marg_h = [0.1, 0.03] ;
    marg_w = [0.1, 0.03] ;
end

% initialize array to store summed force/torque
N_pts = length(t) ;
F_mat = nan(N_pts, 3) ;
T_mat = nan(N_pts, 3) ;

% forces
hatVecForce = eye(3) ;
if plotFlag
    ax_force = gobjects(3) ;
    labelsForce = {'F_{fwd} (N)', 'F_{side} (N)', 'F_{vert} (N)'} ;
end
for i = 1:3
    % unit vector to pull out components
    hat_vec = repmat(hatVecForce(i,:),N_pts,1) ;

    % get right and left wing forces
    F_R = dot(qs_force_struct(1).F_tot, hat_vec,2) ;
    F_L = dot(qs_force_struct(2).F_tot, hat_vec,2) ;
    F_sum = F_R + F_L ;
    
    if plotFlag
        % initialize subplot
        ax_force(i) = subtightplot(3,2,2*(i-1) + 1, gap, marg_h, marg_w) ;
        hold on
        
        % plot components and sum
        plot(t, F_R, '-', 'Color', [1, 0, 0, 0.5],'LineWidth',0.75)
        plot(t, F_L, '-', 'Color', [0, 0, 1, 0.5],'LineWidth',0.75)
        plot(t, F_sum,'-', 'Color', [0, 0, 0, 0.5],'LineWidth',0.75)
        
        % axis properties
        ylabel(labelsForce{i})
        if (i ~= 3)
            ax_force(i).XAxis.Visible = 'off' ;
        else
            xlabel('Time (s)')
        end
        axis tight
        box off
    end
    % store force
    F_mat(:,i) = smooth(F_sum,11,'rloess') ;
end

% torques
% hatVecTorque = [qs_force_struct(1).roll_hat' ; ...
%     qs_force_struct(1).pitch_hat' ; ...
%     qs_force_struct(1).yaw_hat' ] ;
hatVecTorque = eye(3) ;
%hatVecTorque(2,:) = -1.*hatVecTorque(2,:) ; 
if plotFlag
    ax_torque = gobjects(3) ;
    labelsTorque = {'T_{roll} (Nm)', 'T_{pitch} (Nm)', 'T_{yaw} (Nm)'} ;
end
for j = 1:3
    % unit vector to pull out components
    hat_vec = repmat(hatVecTorque(j,:),N_pts,1) ;
    
    % get right and left wing forces
    T_R = dot(qs_force_struct(1).T_tot, hat_vec,2) ;
    T_L = dot(qs_force_struct(2).T_tot, hat_vec,2) ;
    T_sum = T_R + T_L ;
    
    if plotFlag
        % initialize subplot
        ax_torque(j) = subtightplot(3,2,2*j, gap, marg_h, marg_w) ;
        hold on
        % plot components and sum
        plot(t, T_R, '-', 'Color', [1, 0, 0, 0.5],'LineWidth',0.75)
        plot(t, T_L, '-', 'Color', [0, 0, 1, 0.5],'LineWidth',0.75)
        plot(t, T_sum,'-', 'Color', [0, 0, 0, 0.5],'LineWidth',0.75)
        
        % axis properties
        ylabel(labelsTorque{j})
        if (j ~= 3)
            ax_torque(j).XAxis.Visible = 'off' ;
        else
            xlabel('Time (s)')
        end
        axis tight
        box off
    end
    % store torque
    T_mat(:,j) = smooth(T_sum,11,'rloess') ;
end

% ----------------------------------
% plot sum on overlapping axes
if plotFlag
    figure ;
    subplot(2,1,1)
    plot(t, F_mat,'-')
    axis tight
    legend(labelsForce)
    
    subplot(2,1,2)
    plot(t, T_mat,'-')
    axis tight
    legend(labelsTorque)
    
%     % ---------------------------------
%     % compare to measured acceleration
%     if (0)
%         bodyPRY = bodyYPR(:, [2,3,1]) ;
%         [angleVel, angleAccel] = diffBodyEulerAngles(bodyPRY, dt) ;
%         [bodyCM_smooth, bodyVel, bodyAccel] = smoothBodyCM(bodyCM) ;
%         [pitchSmooth, yawSmooth, rollSmooth] = smoothBodyAngles(data) ;
%         pitchSmooth = DEG2RAD.*pitchSmooth ;
%         yawSmooth = DEG2RAD.*yawSmooth ;
%         rollSmooth = DEG2RAD.*rollSmooth ;
%         
%         
%     end
end


end