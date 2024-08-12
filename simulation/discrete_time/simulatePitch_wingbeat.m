% -------------------------------------------------------------------------
% wingbeat time scale simulation of longitudinal flight
% -------------------------------------------------------------------------
function [s_mat, delta_phi_f_all] = ...
    simulatePitch_wingbeat(controllerTerms, t_wb, s_0, thetaB0, params,...
        controlFlag, pertFlag)
% -------------------------
%% read out params
% wing kinematic params
wing_kin_params = read_wing_kin_params(params) ;
omega = wing_kin_params(1) ;
wbf = omega/(2*pi) ;
dt = 1/wbf ;

% morphology
mass = params.body_mass ;           % body mass (kg)
Iyy = params.Iyy ;                  % pitch moment of inertia
C_friction = params.C_friction ;    % rotational drag
g = params.g ;                      % gravitational acceleration (m/s^2)

% general wing kin parameters
wing_kin_params = read_wing_kin_params(params) ;

% body parameters
body_params = read_body_params(params) ; 

% time range for evaluating force/torque in seconds
%  (corresponds to one wingbeat)
t_range = linspace(0, dt, 50)' ; 

% initialize storage for all delta phi front variables
delta_phi_f_all = zeros(size(t_wb)) ; 
% -----------------------------------------
%% controller terms
if controlFlag
    K_i = controllerTerms(1) ;
    K_p = controllerTerms(2) ;
    % NB: time delay should be in wingbeats now!
    deltaT = controllerTerms(3) ;
    if (deltaT < 0.02)
        deltaT = round(deltaT * wbf) ;
    end
else
    K_i = 0 ;
    K_p = 0 ;
    deltaT = 0 ;
end
% ------------------------------------------
%% perturbation info
if pertFlag
    pulseStart = params.pulseStart ;
    pulseEnd = params.pulseEnd ;
    pulseStrength = params.pulseStrength ;
    
    if ~isinteger(pulseEnd)
        pulseEnd = 1 ; 
        % keyboard
    end
end
% pin angle
pinAngle = pi/4 ; % radians. assumed pitch angle of magnetic pin wrt body

% -----------------------------------------------------------
%% iterate over wingbeats
% initialize array to store state information
s_mat = nan(length(t_wb), length(s_0)) ;
s_mat(1:(deltaT+1),:) = repmat(s_0', (deltaT+1), 1) ;

% loop over wingbeats
for ind = (deltaT + 1):(length(t_wb)-1)
    % tic
    % --------------------------------
    % get current body kinematics
    x_dot = s_mat(ind,2) ;
    z_dot = s_mat(ind,4) ; 
    thetaB = s_mat(ind,5) ; 
    thetaB_dot = s_mat(ind, 6) ;
    
    % vectors of velocities
    velBodyBody = [x_dot, 0, z_dot] ;
    bodyYPR_dot = [0, thetaB_dot, 0] ;  %[0, -1*thetaB_dot, 0] ; 
     
    % ---------------------------------
    % calculate PI controller output
    thetaB_lag = s_mat(ind-deltaT, 5) ; 
    thetaB_dot_lag = s_mat(ind-deltaT, 6) ; 
    delta_phi_front = K_i*(thetaB_lag - thetaB0) + K_p*thetaB_dot_lag ; 
    
    % make sure delta_phi_front doesn't grow too large
    delta_phi_front = sign(delta_phi_front)*min([abs(delta_phi_front), ...
        35*(pi/180)]) ; 
    
    % store delta phi front
    delta_phi_f_all(ind+1) = delta_phi_front ; 
    % ----------------------------------------
    % update wing kin params
    wing_kin_params(end) = delta_phi_front ; 
    
    % -------------------------------------------
    % calculate forces and torques at given step
    [F_tot, T_tot] = fakeWingKinForceAndTorque_optim(wing_kin_params, ...
        body_params, t_range, velBodyBody, bodyYPR_dot, params) ;
    
    % get wingbeat mean of forces/torques
    F_mean = mean(F_tot) ; 
    T_mean = mean(T_tot) ; 
    
    % -----------------------------------------------
    %% update equations of motion
    % get the ds/dt vector for current time step (aka wingbeat) 
    ds = zeros(1,6) ; 
    % fwd/back motion
    ds(1) = s_mat(ind, 2) ;
    ds(2) = (1/mass)*F_mean(1) - g*sin(thetaB) + thetaB_dot*z_dot ; % +

    % up/down motion
    ds(3) = s_mat(ind, 4) ; 
    ds(4) = (1/mass)*F_mean(3) - g*cos(thetaB) - thetaB_dot*x_dot ; % -

    % pitching motion
    ds(5) = s_mat(ind, 6) ; 
    ds(6) = (1/Iyy)*(-T_mean(2) - C_friction*thetaB_dot) ;  % -
    
    % add perturbation?
    if pertFlag && (t_wb(ind) > pulseStart) && (t_wb(ind) <= pulseEnd)
        ds(6) = ds(6) + (pulseStrength)*cos(thetaB - pinAngle) ;
    end
    
    % -------------------------------------------------------
    % get prediction for next time step based on current one
    s_mat(ind+1,:) = s_mat(ind,:) + dt.*ds ;
    
    % toc
end

end