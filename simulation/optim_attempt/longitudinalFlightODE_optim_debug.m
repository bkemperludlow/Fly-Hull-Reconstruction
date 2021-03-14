% -------------------------------------------------------------------------
% function to set up longitudinal (i.e. up/down, fwd/back, and
% pitch) flight ordinary differential equation
%
% *** wingbeat averaged?
%   -> doesn't seem to work with ODE solver -- would need to just
%   explicitly solve
%
% INPUTS:
%   - t: time
%   - s: state vector (looks like: s = [x xdot z zdot theta thetadot])
%   - Z: time lagged state vector if solving with dde23. otherwise not used
%   - params: morphological and kinematic parameters
%   - ***Need to add stuff for pert 
%
% OUTPUTS:
%   - ds: change in state vector w.r.t. time
% -------------------------------------------------------------------------
function [ds, delta_phi_f] = longitudinalFlightODE_optim_debug(t, s, Z, ...
    params, thetaB0, pertFlag)
%----------------------------
%% initialize parameters
global delta_phi_f_all
% % body pitch angle set point
% thetaB0 = params.beta_0 ;  % radians
pinAngle = pi/4 ; % radians. assumed pitch angle of magnetic pin wrt body

% controller params
K_i = params.K_i ;                  % integral gain (unitless)
K_p = params.K_p ;                  % proportional gain (seconds)

% morphology
mass = params.body_mass ;           % body mass (kg)
Iyy = params.Iyy ;                  % pitch moment of inertia
C_friction = params.C_friction ;    % rotational drag 
g = params.g ;                      % gravitational acceleration (m/s^2)
% --------------------------------------------
%% read out current and lagged state vectors
% position vars
% x = s(1) ; % fwd/back cm position
% z = s(3) ; % up/down cm position
x_dot = s(2) ; % fwd/back velocity
z_dot = s(4) ; % up/down velocity

% angle vars
thetaB = s(5) ;         % body pitch angle
thetaB_dot = s(6) ;     % body pitch velocity

% switch what we read out depending on whether or not we're using time
% delayed solver
if ~isempty(Z)
    thetaB_lag = Z(5,:) ;     % LAGGED body pitch angle
    thetaB_dot_lag = Z(6,:) ; % LAGGED body pitch velocity
else
    % if we're not solving time-delayed equation, eliminate controller
    % contribution
    thetaB_lag = 0 ; 
    thetaB_dot_lag = 0 ; 
    
    % set controller coeffs to 0
    K_i = 0 ; 
    K_p = 0 ; 
end

% -----------------------------------------------------
%% get current wing and body kinematic params
% change in phi front due to controller
delta_phi_f = K_i.*(thetaB_lag - thetaB0) + K_p.*(thetaB_dot_lag) ;

delta_phi_f_all = [delta_phi_f_all ; delta_phi_f] ; 
% get wing kinematics parameter vector
wing_kin_params = read_wing_kin_params(params, delta_phi_f) ;

% since it's longitudinal flight, left and right wings should be the same
wing_kin_params = repmat(wing_kin_params,1,2) ; 

% also read out body kinematic params
body_params = read_body_params(params) ; 

% ------------------------------------
%% compile body motion into matrices
velBodyBody = [x_dot, 0, z_dot] ;
bodyYPR     = [0, thetaB, 0] ; 
bodyYPR_dot = [0, thetaB_dot, 0] ;  %[0, -1*thetaB_dot, 0] ; 

% -----------------------------------------------------------------
%% calculate forces and torques
% tic
[F_tot, T_tot] = fakeWingKinForceAndTorque_optim(wing_kin_params, ...
    body_params, t, velBodyBody, bodyYPR_dot, params) ;
% toc
% -------------------------------------------------------------------
%% use aerodynamic forces to update equations of motion
%state vector looks like s = [x, x_dot, z, z_dot, theta, theta_dot]
ds = zeros(6,1) ;

% fwd/back motion
ds(1) = s(2) ;
ds(2) = (1/mass)*F_tot(1) - g*sin(thetaB) + thetaB_dot*z_dot ; % +

% up/down motion
ds(3) = s(4) ; 
ds(4) = (1/mass)*F_tot(3) - g*cos(thetaB) - thetaB_dot*x_dot ; % -

% pitching motion
ds(5) = s(6) ; 
ds(6) = (1/Iyy)*(-T_tot(2) - C_friction*thetaB_dot) ; 
%ds(6) = -1*(1/Iyy)*T_tot(2) - C_friction*thetaB_dot ; 
% ------------------------
%% add perturbation?
if pertFlag && (t > params.pulseStart) && (t < params.pulseEnd)
    ds(6) = ds(6) + (params.pulseStrength)*cos(thetaB - pinAngle) ; 
end

% if t > 0.09
%    keyboard 
% end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------
% %% controller equation(s)
% % implements PI control
% function delta_phi_f = ...
%     pitch_controller_func(thetaB_lag, thetaB_dot_lag, thetaB0, K_i, K_p)
% 
% % change in fwd stroke amplitude is determined by body pitch angle and vel
% delta_phi_f = K_i.*(thetaB_lag - thetaB0) + K_p.*(thetaB_dot_lag) ;
% delta_phi_f = nanmean(delta_phi_f) ; 
% 
% if abs(delta_phi_f) > 0.8727
%    delta_phi_f = sign(delta_phi_f)*0.8727 ;
% end
% 
% end

% -----------------------------------------------------------------
%% read out data from params to a wing_kin_params vector
function wing_kin_params = read_wing_kin_params(params, delta_phi_f)
% -----------------------------------------

% list of qs_params wing kin fields
field_list = {'omega', 'phi_f', 'phi_b', 'K', 'theta_0', 'theta_m', ...
    'del_theta','psi_0', 'psi_m', 'del_psi', 'C'} ; 

% factors to convert to radians, if necessary
scale_vec = ones(length(field_list),1) ;  
scale_vec([2,3,5,6,8,9]) = (pi/180).*scale_vec([2,3,5,6,8,9]) ;

% ------------------------------------------
% initialize wing kin params vector
wing_kin_params = zeros(length(field_list)+1,1) ; 

% loop through fields and fill vec
for k = 1:length(field_list)
    wing_kin_params(k) = scale_vec(k)*params.(field_list{k}) ; 
end

% -------------------------------
% add delta phi front
wing_kin_params(length(field_list)+1) = delta_phi_f ; 

end

% -----------------------------------------------------------------
%% read out data from params to a body_params vector
function body_params = read_body_params(params)
% list of qs_params body fields
field_list = {'span', 'r_hinge', 'thorax_width'} ; 

% initialize wing kin params vector
body_params = zeros(length(field_list),1) ; 

% loop through fields and fill vec
for k = 1:length(field_list)
    body_params(k) = params.(field_list{k}) ; 
end

end

% % ------------------------------------------------------------------
% %% functions to get wing Euler angles from parameterized equations
% % ----------------
% % stroke angle
% % ----------------
% function [phi, phi_dot] = phi_func(t, omega, phi_f, phi_b, K)
% % convert forward and backward stroke limits to midpoint and amp
% phi_0 = (phi_f + phi_b)/2 ;
% phi_m = (phi_b - phi_f)/2 ; 
% 
% % get stroke angle and derivative of stroke angle w.r.t. time
% phi = phi_0 + phi_m.*(asin(K.*sin(omega.*t)))/asin(K) ;
% phi_dot = (phi_m*K*omega.*cos(omega.*t))./sqrt(1 - K^2.*(sin(omega*t)).^2) ; 
% 
% end
% 
% % ----------------
% % deviation angle
% % ----------------
% function [theta, theta_dot] = ...
%     theta_func(t, omega, theta_0, theta_m, del_theta)
% 
% % get deviation angle and derivative of deviation angle w.r.t. time
% theta = theta_0 + theta_m.*cos(2*omega.*t + del_theta) ; 
% theta_dot = -2*theta_m*omega.*sin(2*omega.*t + del_theta) ; 
% 
% end
% 
% % ----------------
% % rotation angle
% % ----------------
% function [psi, psi_dot] = psi_func(t, omega, psi_0, psi_m, del_psi, C)
% 
% % get deviation angle and derivative of deviation angle w.r.t. time
% psi = psi_0 + psi_m.*tanh(C.*sin(omega.*t + del_psi))./tanh(C) ; 
% psi_dot = psi_m*C*omega.*cos(omega.*t + ...
%     del_psi).*(sech(C.*sin(omega*t + del_psi))).^2 ;
% 
% end
