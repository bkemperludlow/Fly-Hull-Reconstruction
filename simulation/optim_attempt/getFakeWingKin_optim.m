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
    getFakeWingKin_optim(wing_kin_params, wingSide, t)
% ----------------------------

% structure of indices for params vecs
% ind_struct = defineParamIndices() ;

% --------------------------------------------
%% get current wing kin params
omega = wing_kin_params(1) ; % wbf (Hz)

phi_f = wing_kin_params(2) ; % stroke angle
phi_b = wing_kin_params(3) ;
K = wing_kin_params(4) ;

theta_0 = wing_kin_params(5) ; % deviation angle
theta_m = wing_kin_params(6) ;
del_theta = wing_kin_params(7) ;

psi_0 = wing_kin_params(8) ; % rotation angle
psi_m = wing_kin_params(9) ;
del_psi = wing_kin_params(10) ;
C = wing_kin_params(11) ;

delta_phi_f = wing_kin_params(12) ;
% if length(wing_kin_params) < 12
%     delta_phi_f = 0 ;
% else
%     delta_phi_f = wing_kin_params(ind_struct.delta_phi_f) ;
% end
% --------------------------------------------
%% use params to get wing kin
% -------------------
% stroke angle
% -------------------
phi_f_curr = phi_f + delta_phi_f ;
% convert forward and backward stroke limits to midpoint and amp
phi_0 = (phi_f_curr + phi_b)/2 ;
phi_m = (phi_b - phi_f_curr)/2 ;

% get stroke angle and derivative of stroke angle w.r.t. time
phi = phi_0 + phi_m.*(asin(K.*sin(omega.*t)))/asin(K) ;
phi_dot = (phi_m*K*omega.*cos(omega.*t))./sqrt(1 - K^2.*(sin(omega*t)).^2) ;

% -------------------
% deviation angle
% -------------------
theta = theta_0 + theta_m.*cos(2*omega.*t + del_theta) ;
theta_dot = -2*theta_m*omega.*sin(2*omega.*t + del_theta) ;

% -------------------
% rotation angle
% -------------------
psi = psi_0 + psi_m.*tanh(C.*sin(omega.*t + del_psi))./tanh(C) ;
psi_dot = psi_m*C*omega.*cos(omega.*t + ...
    del_psi).*(sech(C.*sin(omega*t + del_psi))).^2 ;


% --------------------------------------------
%% adjust based on wing side
if contains(wingSide,'R','IgnoreCase',true)
    % sign flip of stroke angle
    phi = -1.*phi ;
    phi_dot = -1.*phi_dot ;
    
    % rotation angle flip -- would normally do this on left wing, but we're
    % giving psi a NEGATIVE phase offset, which means it needs to happen on
    % the right wing
    psi_dot = -1.*psi_dot ;
    psi = pi - psi ;
end

% --------------------------------------------------------
%% put euler angles and derivatives into output matrices
wingAngleMat = [phi, theta, psi] ;
wingVelMat = [phi_dot, theta_dot, psi_dot] ;


end


