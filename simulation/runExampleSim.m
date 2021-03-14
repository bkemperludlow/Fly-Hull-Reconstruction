% -------------------------------------------------------------------------
% script to run a simulation for uncontrolled flight and plot results
% -------------------------------------------------------------------------
%-----------------------
%% params
% input flags
pertFlag = true ;
controlFlag = true ;
plotFlag = true ;

% ----------------------------------------------------------
% general params struct
params = defineQuasiSteadyParams ;
omega = params.omega ; 

% ----------------------------------------
% initial state guess
thetaB0     =   pi/4 ; % rad
x_0         =   0 ; % m
xdot_0      =   0 ; % m/s
z_0         =   0 ; % m
zdot_0      =   0 ; % m/s
theta_0     =   pi/4 ; % rad
thetadot_0  =   0 ;  % rad/s
s_0 = [x_0; xdot_0; z_0; zdot_0; theta_0 ; thetadot_0] ;

% ----------------------------------
% type of delay function to use
% 'constant' = just deltaT (updates at each time point)
% 'average'  = average over previous wingbeat, also different at every time point
% 'function' = use "myControllerDelay.m" to get delay. should give same
%              value for multiple successive time points
delayFunType = 'function' ;  % 'constant' | 'average' | 'function'

% --------------------------------
% limit time of simulation run?
timeLimSec = 80 ; % nan ; %80 ;

% -------------------------------------------------------
% controller terms
% ---------------------
% OLD
K_i = 0.537556 ; % 0.437556 ;
K_p = 0.005728 ; % 0.009728 ; 
deltaT = 4.1739e-3; % 4.3478e-3 ; %0.00430 ;
% ---------------------

% ---------------------------------
% OPTIMIZED (with delay function)
% K_i = 0.573607 ;
% K_p = 0.005218 ;
% deltaT = 0.004574 ;
% ----------------------------------

% -------------------------------------------------------
% OPTIMIZED (with time delay averaged over wingbeat)
% K_i = 0.633928 ;
% K_p = 0.006656 ;
% deltaT = 0.004491 ;
% -------------------------------------------------------

% % -----------------------------------------
% % OPTIMIZED (10/10/2020)
% K_i = 0.408409 ; % 0.408409 ;
% K_p = 0.009572 ;
% deltaT = 0.004449 ;
% % -----------------------------------------

% -----------------------------------------
% OPTIMIZED (10/12/2020)
% K_i = 0.859873 ;
% K_p = 0.006050 ;
% deltaT = 0.004373 ;
% -----------------------------------------

% ---------------------
% TEST
% K_i = 0 ; 
% K_p = 0.01218 ;
% deltaT = 0.008574 ;
% ---------------------

% ---------------------------------------------
% OPTIMIZED (10/14/20, with function delay)
% K_i = 0.379619 ;
% K_p = 0.018663/3 ;
% deltaT = 0.004766 ; 
% ---------------------------------------------

controllerTerms = [K_i, K_p, deltaT] ;

% -------------------------------------------------------------
% time interval
% start point
t_int1 = -0.25 ;  
% t_int1 = t_int1 + (-1*rem(omega*t_int1, 2*pi) + pi/2)/omega; % ensure we start with sin(t_int1) = 1 
t_int1 = t_int1 + (-1*rem(omega*t_int1, 2*pi) - pi/2 - pi/8)/omega; % ensure we start with sin(t_int1) = 1 
% end point
t_int2 = 0.5 ; 

% interval(s)
t_int = [t_int1, t_int2] ;  % t_int = [-0.15, tf + 0.050] ; 
% t_int = (1e-3).*[1.1, 200] ; % (1e-3).*[1.1, 1000] ;  %(1e-3).*[-200, 500] ;

% plot params
figPosition =  [202, 389, 1046, 674] ;

% ----------------------------------------------------
%% run simulation
sol = simulatePitch_optim(controllerTerms, t_int, s_0, thetaB0, ...
    params, controlFlag, pertFlag, delayFunType, timeLimSec) ;
% sol = simulatePitch_optim_debug(controllerTerms, t_int, s_0, thetaB0, ...
%     params, controlFlag, pertFlag) ;
% ----------------------------------------------------
%% read out data from solution struct
% ----------------------------------
% misc params
omega = params.omega ;
time_label = 'Time (ms)' ;

% ----------------------------------
% read out data from solution
dt = 0.125e-3 ; % seconds
t_eval = (min(sol.x) : dt : max(sol.x)) ;

if isfield(sol,'yp')
    [sint, dsint] = deval(sol, t_eval) ; %evaluate solution at tspan
else
    sint = deval(sol, t_eval) ; %evaluate solution at tspan
end

thetaB = sint(5,:) ; %radians
thetaB_dot = sint(6,:) ; % radians/sec
x = sint(1,:) ; % meters
z = sint(3,:) ; % meters
x_dot = sint(2,:) ; % m/s
z_dot = sint(4,:) ; % m/s

% ------------------------------------------------
%% change units/frame
% get lab frame x and z
[x_lab, x_dot_lab, z_lab, z_dot_lab] = ...
    convertSimCoordsToLab(x_dot, z_dot, thetaB, dt) ;

% change units
x_lab = 1000.*x_lab ; % m -> mm
z_lab = 1000.*z_lab ; % m -> mm
%     x_dot_lab = 1000.*x_dot_lab ; % m/s -> mm/s
%     z_dot_lab = 1000.*z_dot_lab ; % m/s -> mm/s
thetaB = (180/pi).*thetaB ; % rad -> deg
thetaB_dot = (180/pi).*thetaB_dot ; % rad/s -> deg/s
t_eval = 1000.*t_eval ;  % s -> ms

% ------------------------------------------------
%% filter body kinematics to see general trend
d1 = designfilt('lowpassiir','FilterOrder',3,'SampleRate',(1/dt), ...
    'HalfPowerFrequency',omega/(8*pi),'DesignMethod','butter');
thetaB_filt = filtfilt(d1, thetaB) ;
thetaB_dot_filt = filtfilt(d1, thetaB_dot) ;

if exist('dsint', 'var')
   pitch_torque = dsint(6,:) ; 
   pitch_torque_filt = filtfilt(d1, pitch_torque) ; 
end
    

% ----------------------------------------------------
%% plot results?
if plotFlag
    h_main = figure('OuterPosition', figPosition) ;
    % ---------------------------------------------
    % real space tracjectory
    
    sub_samp = 10 ; % sub sample results to keep plot from getting crowded
    ind = 1:sub_samp:length(t_eval) ;
    
    % initialize figure
    ax_coords = subplot(2,2,1) ;
    hold on
    
    % plot locations colored by pitch angle
    scatter(x_lab(ind), z_lab(ind), [], thetaB(ind), 'o', 'filled')
    phasemap
    
    % axis properties
    axis tight
    axis equal
    xlabel('x (mm)')
    ylabel('z (mm)')
    
    % ---------------------------------------------
    % real space velocity
    ax_vel = subplot(2,2,3) ;
    hold on
    
    hvelx = plot(t_eval, 100.*x_dot_lab, 'k-') ;
    hvelz = plot(t_eval, 100.*z_dot_lab, 'r-') ;
    
    % axis properties
    axis tight
    xlabel(time_label)
    ylabel('velocity (cm/s)')
    
    legend([hvelx, hvelz], {'x','z'}, 'location', 'northwest')
    
    % ---------------------------------------------
    % body pitch vs time
    ax_pitch = subplot(2,2,2) ;
    hold on
    %     plot(t_eval, thetaB - (180/pi)*thetaB0, '-')
    %     plot(t_eval, thetaB_filt - (180/pi)*thetaB0, 'r-','LineWidth',2)
    plot([t_eval(1), t_eval(end)], (180/pi)*thetaB0.*[1, 1], 'k--',...
        'LineWidth',1.25)
    plot(t_eval, thetaB, '-', 'LineWidth',0.5)
    plot(t_eval, thetaB_filt, 'r-','LineWidth',1.5)

    % axis properties
    axis tight
    xlabel(time_label)
    ylabel('Body Pitch Angle (deg)')

    % ---------------------------------------------
    % body pitch velocity vs time
    ax_pitch_vel = subplot(2,2,4) ;
    hold on
    plot(t_eval, thetaB_dot, '-')
    plot(t_eval, thetaB_dot_filt, 'r-','LineWidth',2)
    
    % axis properties
    axis tight
    xlabel(time_label)
    ylabel('Body Pitch Vel (deg/s)')
    
end