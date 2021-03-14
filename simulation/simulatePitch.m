% -------------------------------------------------------------------------
% updated function to run longitudinal flight simulations to test control
% parameters
%
% INPUTS:
%   - controllerTerms: 3-element vector specifying PI controller quantities
%       of the form [K_i, K_p, deltaT]
%   - t_int: interval over which to run de solver, specified as 2-element
%       vector. In seconds
%   - s_0: initial conditions state vector (also 'history' for dde). of the
%       form s_0 = [x, x_dot, z, z_dot, thetaB, thetaB_dot]
%   - thetaB0: controller angle set point
%   - params: quasi steady param structure, as defined by
%       defineQuasiSteadyParams.m
%   - controlFlag: boolean specifying whether or not to implement control
%   - pertFlag: boolean indicating whether or not to apply perturbation.
%       pert parameters are stored in params structure.
%   - plotFlag: boolean to plot results or not
%
% OUTPUT:
%   - sol: solution output structure from DE solver. can use deval.m to get
%       values of solution at arbitrary time points in range
%
% -------------------------------------------------------------------------
function sol = simulatePitch(controllerTerms, t_int, s_0, thetaB0, params, ...
    controlFlag, pertFlag, plotFlag)
%-----------------------
%% inputs
if ~exist('controllerTerms','var') || isempty(controllerTerms)
%     % controller base values (from JEB paper)
%     K_i = 0.3 ; %0.3 ; %.3 ;
%     K_p = 0.007 ; %.007 ; %seconds
%     deltaT = 0.006 ; %seconds %0.006
%     controllerTerms = [K_i, K_p, deltaT] ;
    % controller base values from fmincon
    K_i = 0.437556 ;
    K_p = 0.005728 ;
    deltaT = 0.002148 ;
    controllerTerms = [K_i, K_p, deltaT] ;
end
if ~exist('t_int','var') || isempty(t_int)
    t_int = (1e-3).*[-99, 200] ;  %(1e-3).*[-200, 500] ; 
end
if ~exist('s_0','var') || isempty(s_0)
    % initial conditions
    x_0         =   0 ; % m
    xdot_0      =   0 ; % m/s
    z_0         =   0 ; % m
    zdot_0      =   0 ; % m/s
    theta_0     =   pi/4 ; % rad
    thetadot_0  =   0 ;  % rad/s
    s_0 = [x_0; xdot_0; z_0; zdot_0; theta_0; thetadot_0] ;
end
if ~exist('thetaB0','var') || isempty(thetaB0)
    thetaB0 = pi/4  ; % radians 
end
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams ;
end
if ~exist('controlFlag','var') || isempty(controlFlag)
    controlFlag = true ; % true ; 
end
if ~exist('pertFlag','var') || isempty(pertFlag)
    pertFlag = false ; % true
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true ;
end

% ----------------------------------------------
%% add controller terms to params struct
if controlFlag
    params.K_i = controllerTerms(1) ;
    params.K_p = controllerTerms(2) ;
    params.deltaT = controllerTerms(3) ;
else
    params.K_i = 0 ;
    params.K_p = 0 ;
    params.deltaT = 0 ;
end

% ----------------------------------------------
%% run solver
if controlFlag && (params.deltaT > 0)
    func = @(t,s,Z) longitudinalFlightODE(t,s,Z,params,thetaB0,pertFlag) ; 
    sol = dde23(func, params.deltaT, s_0, t_int) ; % t_int
else
    % NB: want to use ode23 here (despite it being inferior to ode45 for
    % our purposes) for consistency with conditions for delayed vs
    % non-delayed cases
    
    func = @(t,s) longitudinalFlightODE(t,s,[],params, thetaB0, pertFlag) ;
    %     sol = ode45(@(t,s) longitudinalFlightODE(t,s,[],params), t_int, s_0) ;
    sol = ode23(func, t_int, s_0) ;
end

% ------------------------------------------------
%% plot results?
if plotFlag
    % ----------------------------------
    % misc params
    omega = params.omega ;
    time_label = 'Time (ms)' ; 
    
    % ----------------------------------
    % read out data from solution
    dt = 1e-4 ; % seconds
    t_eval = (min(sol.x) : dt : max(sol.x)) ;
    
    sint = deval(sol, t_eval) ; %evaluate solution at tspan
    thetaB = sint(5,:) ; %radians
    thetaB_dot = sint(6,:) ; % radians/sec
    x = sint(1,:) ; % meters
    z = sint(3,:) ; % meters
    x_dot = sint(2,:) ; % m/s
    z_dot = sint(4,:) ; % m/s
    
    % ------------------------------------------------
    %% change units/frame
    % get lab frame x and z
    x_lab = zeros(size(x)) ;
    z_lab = zeros(size(z)) ;
    x_dot_lab = zeros(size(x_dot)) ;
    z_dot_lab = zeros(size(z_dot)) ;
    
    for k = 1:length(x)
        % current rotation matrix
        rotM = eulerRotationMatrix(0, thetaB(k),0) ;
        
        % rotate position vector
        pos_vec = [x(k) ; 0 ; z(k)] ;
        pos_vec_lab = rotM'*pos_vec ;
        x_lab(k) = pos_vec_lab(1) ;
        z_lab(k) = pos_vec_lab(3) ;
        
        % rotate velocity vector
        vel_vec = [x_dot(k) ; 0 ; z_dot(k)] ; 
        vel_vec_lab = rotM'*vel_vec ; 
        x_dot_lab(k) = vel_vec_lab(1) ; 
        z_dot_lab(k) = vel_vec_lab(3) ; 
    end
    
    % change units
    % x = 1000.*x ; % m -> mm
    % z = 1000.*z ; % m -> mm
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
        'HalfPowerFrequency',omega/(4*pi),'DesignMethod','butter');
    thetaB_filt = filtfilt(d1, thetaB) ;
    thetaB_dot_filt = filtfilt(d1, thetaB_dot) ;
    
    
    % ---------------------------------------------
    % REAL SPACE TRAJECTORY
    % ---------------------------------------------
    % sub sample results to keep plot from getting crowded
    sub_samp = 10 ;
    ind = 1:sub_samp:length(t_eval) ;
    
    % initialize figure
    h1 = figure ;
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
    % BODY PITCH VS TIME
    % ---------------------------------------------
    h2 = figure ;
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
    % BODY PITCH VEL VS TIME
    % ---------------------------------------------
    h3 = figure ;
    hold on
    plot(t_eval, thetaB_dot, '-')
    plot(t_eval, thetaB_dot_filt, 'r-','LineWidth',2)
    
    % axis properties
    axis tight
    xlabel(time_label)
    ylabel('Body Pitch Vel (deg/s)')
    
    
%     % ---------------------------------------------
%     % TRANSLATIONAL VEL VS TIME
%     % ---------------------------------------------
%     h4 = figure ;
%     yyaxis left
%     plot(t_eval, z_dot)
%     yyaxis right
%     plot(t_eval, thetaB)
%     yyaxis left
%     hold on
%     plot(t_eval, x_dot)
     
    % axis properties
    %     % ----------------------------------------------
    %     % LATE-TERM PITCH VS PITCH VEL
    %     % ----------------------------------------------
    %     t_late_ind = t_eval > 50 ;
    %     h4 = figure ;
    %     plot(thetaB(t_late_ind) - thetaB0, thetaB_dot(t_late_ind))
    %     xlabel('Body Pitch Angle (deg)')
    %     ylabel('Body Pitch Vel. (deg/s)')
end
end