% -------------------------------------------------------------------------
% function to get wing euler angles and their derivatives w.r.t. time for a
% a set of parameters (according to the functional forms described in
% (Chang and Wang, PNAS, 2014))
%
% INPUTS:
%   - wing_kin_params: 12x1 vector containing parameter values for fake
%     wing kinenmatics of the form:
%       wing_kin_params = [omega; phi_f; phi_b; K; theta_0; theta_m;...
%                             del_theta; psi_0; psi_m; del_psi; C; ...
%                             delta_phi_front]
%   - wingSide: string ('R' | 'L') for right or left wing, respectively
%   - t: Nx1 vector of time values over which to evaluate the force/torque
%
% OUTPUTS:
%   - wingAngleMat: Nx3 matrix of wing Euler angles of the form [phi, ...
%       theta, psi (aka eta)]
%   - wingVelMat: Nx3 matrix of wing Euler angles derivatives w.r.t. time.
%       same form as above
% -------------------------------------------------------------------------
function [wingAngleMat, wingVelMat] = ...
    getFakeWingKin(wing_kin_params, wingSide, t, plotFlag)
% ----------------------------
%% inputs and params
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = false ;
end

% structure of indices for params vecs
ind_struct = defineParamIndices() ;

% --------------------------------------------
%% get current wing kin params
omega = wing_kin_params(ind_struct.omega) ; % wbf (Hz)

phi_f = wing_kin_params(ind_struct.phi_f) ; % stroke angle
phi_b = wing_kin_params(ind_struct.phi_b) ;
K = wing_kin_params(ind_struct.K) ;

theta_0 = wing_kin_params(ind_struct.theta_0) ; % deviation angle
theta_m = wing_kin_params(ind_struct.theta_m) ;
del_theta = wing_kin_params(ind_struct.del_theta) ;

psi_0 = wing_kin_params(ind_struct.psi_0) ; % deviation angle
psi_m = wing_kin_params(ind_struct.psi_m) ;
del_psi = wing_kin_params(ind_struct.del_psi) ;
C = wing_kin_params(ind_struct.C) ;

if length(wing_kin_params) < 12
    delta_phi_f = 0 ;
else
    delta_phi_f = wing_kin_params(ind_struct.delta_phi_f) ;
end
% --------------------------------------------
%% use params to get wing kin
% stroke angle
phi_f_curr = phi_f + delta_phi_f ; 
[phi, phi_dot] = phi_func(t, omega, phi_f_curr, phi_b, K) ;

% deviation angle
[theta, theta_dot] = theta_func(t, omega, theta_0, theta_m, del_theta) ;

% rotation angle
[psi, psi_dot] =  psi_func(t, omega, psi_0, psi_m, del_psi, C) ;

% --------------------------------------------
%% adjust based on wing side
if contains(wingSide,'R','IgnoreCase',true) && (del_psi < 0)
    % sign flip of stroke angle
    phi = -1.*phi ;
    phi_dot = -1.*phi_dot ;
    
    % rotation angle flip
    psi_dot = -1.*psi_dot ;
    psi = pi - psi ;
elseif contains(wingSide,'R','IgnoreCase',true) && (del_psi > 0)
    % sign flip of stroke angle
    phi = -1.*phi ;
    phi_dot = -1.*phi_dot ;
elseif contains(wingSide,'L','IgnoreCase',true) && (del_psi > 0)
    % rotation angle flip
    psi_dot = -1.*psi_dot ;
    psi = pi - psi ;
else
    fprintf('Invalid wing side: %s \n', wingSide)
    keyboard
end

% --------------------------------------------------------
%% put euler angles and derivatives into output matrices
% make sure output is Nx3 (rather than 3xN)
if size(phi,1) < size(phi,2)
    % in this case, phi is a row vector, and need to flip
    wingAngleMat = [phi', theta', psi'] ;
    wingVelMat = [phi_dot', theta_dot', psi_dot'] ;
else
    % in this case, everything is already a column vector (as it should be)
    wingAngleMat = [phi, theta, psi] ;
    wingVelMat = [phi_dot, theta_dot, psi_dot] ;
end

% --------------------------------------------------------
%% plot results?
if plotFlag
    h_main = figure ;
    angleLabels = {'Stroke', 'Deviation', 'Rotation'} ;
    cc = 1 ;
    for dim = 1:3
        % ----------------
        % plot angle
        % ----------------
        subplot(3,2,cc)
        plot(t, (180/pi).*wingAngleMat(:,dim))
        
        % axis
        xlabel('Time (sec)')
        ylabel([angleLabels{dim} ' Angle (deg)'])
        axis tight
        
        % increment counter 
        cc = cc + 1 ; 
        
        % ----------------
        % plot velocity
        % ----------------
        subplot(3,2,cc)
        plot(t, (180/pi).*wingVelMat(:,dim))
        
        % axis
        xlabel('Time (sec)')
        ylabel([angleLabels{dim} ' Velocity (deg/s)'])
        axis tight
        
        % increment counter 
        cc = cc + 1 ; 
        
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------
%% functions to get wing Euler angles from parameterized equations
% ----------------
% stroke angle
% ----------------
function [phi, phi_dot] = phi_func(t, omega, phi_f, phi_b, K)
% convert forward and backward stroke limits to midpoint and amp
phi_0 = (phi_f + phi_b)/2 ;
phi_m = (phi_b - phi_f)/2 ;

% get stroke angle and derivative of stroke angle w.r.t. time
phi = phi_0 + phi_m.*(asin(K.*sin(omega.*t)))/asin(K) ;
phi_dot = (phi_m*K*omega.*cos(omega.*t))./sqrt(1 - K^2.*(sin(omega*t)).^2) ;

end

% ----------------
% deviation angle
% ----------------
function [theta, theta_dot] = ...
    theta_func(t, omega, theta_0, theta_m, del_theta)

% get deviation angle and derivative of deviation angle w.r.t. time
theta = theta_0 + theta_m.*cos(2*omega.*t + del_theta) ;
theta_dot = -2*theta_m*omega.*sin(2*omega.*t + del_theta) ;

end

% ----------------
% rotation angle
% ----------------
function [psi, psi_dot] = psi_func(t, omega, psi_0, psi_m, del_psi, C)

% get deviation angle and derivative of deviation angle w.r.t. time
psi = psi_0 + psi_m.*tanh(C.*sin(omega.*t + del_psi))./tanh(C) ;
psi_dot = psi_m*C*omega.*cos(omega.*t + ...
    del_psi).*(sech(C.*sin(omega*t + del_psi))).^2 ;

end

% ----------------------------------------------------
%% define indices for param vecs
function ind_struct = defineParamIndices()
% write out indices in wing params vec
omega_ind = 1 ;

phi_f_ind = 2 ;
phi_b_ind = 3 ;
K_ind = 4 ;

theta_0_ind = 5 ;
theta_m_ind = 6 ;
del_theta_ind = 7 ;

psi_0_ind = 8 ;
psi_m_ind = 9 ;
del_psi_ind = 10 ;
C_ind = 11 ;

delta_phi_f_ind = 12 ; 

% write out indices in body params vec
span_ind = 1 ;
hinge_vec_ind = 2 ;
thorax_width_ind = 3 ;

% add indices to struct
ind_struct = struct() ;
ind_struct.omega = omega_ind ;

ind_struct.phi_f = phi_f_ind ;
ind_struct.phi_b = phi_b_ind ;
ind_struct.K = K_ind ;

ind_struct.theta_0 = theta_0_ind ;
ind_struct.theta_m = theta_m_ind ;
ind_struct.del_theta = del_theta_ind ;

ind_struct.psi_0 = psi_0_ind ;
ind_struct.psi_m = psi_m_ind ;
ind_struct.del_psi = del_psi_ind ;
ind_struct.C = C_ind ;

ind_struct.delta_phi_f = delta_phi_f_ind ;

ind_struct.span = span_ind ;
ind_struct.hinge_vec = hinge_vec_ind ;
ind_struct.thorax_width = thorax_width_ind ;

end