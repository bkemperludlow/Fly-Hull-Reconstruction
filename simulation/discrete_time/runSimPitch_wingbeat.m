% -------------------------------------------------------------------------
% script to test wing-beat averaged simulation of longitudinal flight
% -------------------------------------------------------------------------
%% params
plotFlag = true ;       % plot results?

controlFlag = true ;   % use PI control?
pertFlag = true ;      % apply perturbation

% time range (in wingbeats)
t_wb = (-20:40)' ;

% general params
params = defineQuasiSteadyParams ;
omega = params.omega ;
dt = (2*pi)/omega ;

% --------------------------------------------
%% define initial conditions
thetaB0     =   pi/4 ; % rad
x_0         =   0 ; % m
xdot_0      =   0 ; % m/s
z_0         =   0 ; % m
zdot_0      =   0 ; % m/s
theta_0     =   pi/4 ; % rad
thetadot_0  =   0 ;  % rad/s
s_0 = [x_0; xdot_0; z_0; zdot_0; theta_0 ; thetadot_0] ;

% ----------------------------------------------
%% define controller terms
K_i = 0.908409 ; % 0.408409 ;
K_p = 0.009572 ;
deltaT = 1 ; % wingbeats!

controllerTerms = [K_i, K_p, deltaT] ; 
% ----------------------------------------------
%% update pulse start/stop info
params.pulseStart = 0 ;
params.pulseEnd = 1 ;

% --------------------------------------------------------------
%% run simulation
[s_mat, delta_phi_f] = simulatePitch_wingbeat(controllerTerms, t_wb, s_0, ...
    thetaB0, params, controlFlag, pertFlag) ;

% -------------------------------------
%% plot results?
if plotFlag
    % read out kinematics from state vector matrix
    x_dot = s_mat(:, 2) ;
    z_dot = s_mat(:, 4) ;
    thetaB = s_mat(:, 5) ;
    thetaB_dot = s_mat(:, 6) ;
    
    % convert x,z to lab frame
    [x_lab, x_dot_lab, z_lab, z_dot_lab] = ...
        convertSimCoordsToLab(x_dot, z_dot, thetaB, dt) ;
    
    % ------------------------------------------------------
    %% make plots
    figPosition =  [202, 389, 1046, 674] ;
    h_main = figure('OuterPosition', figPosition) ;
    time_label = 'Time (wb)' ; 
    
    % ---------------------------------------------
    % real space tracjectory
    ax_coords = subplot(2,2,1) ;
    hold on
    
    % plot locations colored by pitch angle
    scatter(1000.*x_lab, 1000.*z_lab, [], thetaB, 'o', 'filled')
    % phasemap
    
    % axis properties
    axis tight
    axis equal
    xlabel('x (mm)')
    ylabel('z (mm)')
    
    % ---------------------------------------------
    % real space velocity
    ax_vel = subplot(2,2,3) ;
    hold on
    
    plot(t_wb, 100.*x_dot_lab, 'k-')
    plot(t_wb, 100.*z_dot_lab, 'r-')
    
    % axis properties
    axis tight
    xlabel(time_label)
    ylabel('velocity (cm/s)')
    
    % ---------------------------------------------
    % body pitch vs time
    ax_pitch = subplot(2,2,2) ;
    hold on
    %     plot(t_eval, thetaB - (180/pi)*thetaB0, '-')
    %     plot(t_eval, thetaB_filt - (180/pi)*thetaB0, 'r-','LineWidth',2)
    plot([t_wb(1), t_wb(end)], (180/pi)*thetaB0.*[1, 1], 'k--',...
        'LineWidth',1.25)
    plot(t_wb, (180/pi).*thetaB, '-', 'LineWidth',1.0)
    %plot(t_eval, thetaB_filt, 'r-','LineWidth',1.5)

    % axis properties
    axis tight
    xlabel(time_label)
    ylabel('Body Pitch Angle (deg)')

    % ---------------------------------------------
    % body pitch velocity vs time
    ax_pitch_vel = subplot(2,2,4) ;
    hold on
    plot(t_wb, (180/pi).*thetaB_dot, '-')
    % plot(t_eval, thetaB_dot_filt, 'r-','LineWidth',2)
    
    % axis properties
    axis tight
    xlabel(time_label)
    ylabel('Body Pitch Vel (deg/s)')
    
    % -----------------------------------------------
    %% also plot delta phi front
    h_phif = figure ; 
    yyaxis left
    plot(t_wb, (180/pi).*delta_phi_f, 'k.-')
    ylabel('\Delta\phi front (deg)')
    
    yyaxis right
    plot(t_wb, (180/pi).*thetaB, '-')
    ylabel('Body Pitch Angle (deg)')
    axis tight
    xlabel('Time (wb)') 
    
    
end
