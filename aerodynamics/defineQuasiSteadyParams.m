% -------------------------------------------------------------------------
% function to generate a structure containing parameters needed for
% calculating aerodynamic forces using the quasi-steady aerodynamic model
% -------------------------------------------------------------------------
function qs_params = defineQuasiSteadyParams(varargin)
% -------------------------------------------
%% parameters upon which other things depend
defaults = {'span', 0.002125 , ...    % meters % was 0.0025
            'r22_S', .313 , ...     % unitless
            'r_hinge', 0.000219 } ;   % meters % was 0.0003 %0.000219

% wing span
span = get_var('span', 'defaults', defaults, varargin{:}); 
% second moment of wing inertia. Sane+Dickinson (2002) gives .4, Cheng et al. (2009) gives .313
r22_S = get_var('r22_S', 'defaults', defaults, varargin{:});
% distance along long body axis from CM to wing hinge
r_hinge = get_var('r_hinge', 'defaults', defaults, varargin{:});

% ------------------------------------
%% physical constants
% gravitational accel
g = 9.8 ; %m/s^2 

% density of air
airDensity = 1.2041 ; %kg/m^3

% ---------------------------------------
%% computational methods
fitType = 'smoothingspline' ; % 'cubicspline' | 'smoothingspline' (for alpha)
wingKinDiffType = 'fit' ; % how to perform differentiation
bodyVelSmoothWin = 51 ; % window size for smoothing body frame velocity

% ------------------------------------
%% wing morphological parameters
% wingbeat frequency in Hz
freq = 225 ; % 230 %Hz %was 250

% wingbeat frequency in rad/s
omega = 2*pi*freq ; 

% % wing span
% span = .0025 ; %meters %.002

% wing chord
chord = .0007 ; %meters %.001

% wing area
S = 2.023e-6 + span*1.748e-3 ; % meters^2, formula from (Fry, Sayaman, Dickinson, 2005)
%S = pi*(span/2)*(chord/2) ; %area of ellipse

% % second moment of wing inertia
% r22_S = .313 ; %Sane+Dickinson (2002) gives .4, Cheng et al. (2009) gives .313

% other wing moments (can derive other wing shape parameters based on 
%  r22_S (see Ellington (1984) and  Whitney & Wood (2010))
r2_S = sqrt(r22_S) ; % non-dimensional radius of 2nd moment of wing area
r1_S = (r2_S/0.929).^(1/0.732) ; % non-dimensional radius of 1st moment of wing area
r12_S = r1_S.^2 ; % squared 1st moment of wing area

% parameters for beta distribution describing wing shape
% NB: c_hat = ( r_hat^(p_beta - 1) * (1 - r_hat)^(q_beta - 1) ) / beta(p_beta, q_beta)
p_beta = r1_S*((r1_S*(1-r1_S))/(r22_S - r12_S) - 1)  ; 
q_beta = (1-r1_S)*((r1_S*(1-r1_S))/(r22_S - r12_S) - 1)  ; 
Beta_pq = beta(p_beta, q_beta) ; 

% integral of area moment used for rotation force: 
%    integral from 0 to 1 of r_hat*c_hat^2 dr_hat
% From wolfram alpha, this integral amounts to :
%    (1/Beta_pq^2)*gamma(2*p_beta)*gamma(2*q_beta - 1) / gamma(2*p_beta + 2*q_beta -1)
c2_r_int = ((1/Beta_pq^2)*gamma(2*p_beta)*gamma(2*q_beta - 1)) / ...
    gamma(2*p_beta + 2*q_beta - 1) ; 

% integral of area moment used for rotation force: 
%    integral from 0 to 1 of c_hat^2 dr_hat
% From wolfram alpha, this integral amounts to :
%    (1/Beta_pq^2)*gamma(2*p_beta-1)*gamma(2*q_beta - 1) / gamma(p_beta + q_beta -1)
c2_int = ((1/Beta_pq^2)*gamma(2*p_beta-1)*gamma(2*q_beta - 1)) / ...
    gamma(2*(p_beta + q_beta -1)) ;

mean_chord = chord/1.2681 ; % using the relation c_hat = chord/mean_chord and finding c_hat maximum
% ------------------------------------
%% body morphological parameteres
% body mass
body_mass = 2.796e-7 + span*4.982e-4 ; % kg, (Fry, Sayaman, Dickinson, 2005) 
%body_mass = 1.1e-6 ; %kg, (Cheng, et al., 2009)

% body length -- should check this
body_length = .0024 ; %meters
%body_width = .0012 ; %meters
thorax_width = 0.000783/2 ; % meters % NB: factor of 1/2 => this is an ellipsoid radius

% moment of inertia (from Cheng (2009)) 
MOI_pitch = .506e-12 ; %N*m*s^2
MOI_roll = .1145e-12 ; %N*m*s^2
MOI_yaw = .4971e-12 ; %N*m*s^2

MOI_rollYaw = -0.191e-12 ; % N*m*s^2

% alternative set for moment of inertia, i think for body frame(?) 
% (from Cheng and Deng, 2011, http://dx.doi.org/10.1109/TRO.2011.2156170)
Ixx = 3.06e-13 ; % N*m*s^2
Iyy = 5.06e-13 ; % N*m*s^2
Izz = 3.06e-13 ; % N*m*s^2
Ixz = -1.91e-13 ; % N*m*s^2

% trying stroke plane moment of inertia from Muijres et al. 2015 
%http://dx.doi.org/10.1242/jeb.114280
f_steady = 189 ; % Hz
Ixx_strokePlane = (0.64*body_mass*g*span)/f_steady^2 ; % N*m*s^2
Iyy_strokePlane = (1.07*body_mass*g*span)/f_steady^2 ; % N*m*s^2
Izz_strokePlane = (0.57*body_mass*g*span)/f_steady^2 ; % N*m*s^2
% resting pitch angle in flight (also stroke plane relative to body angle)
beta_0 = pi/4 ; %radians

% ------------------------------------
%% stroke plane and wing base
%rotates the velocities and such into the stroke plane in the body axis
R_sp = [ cos(pi/2-beta_0),     0, sin(pi/2-beta_0) ; ...
                          0,     1,                  0 ; ...
        (-1)*sin(pi/2-beta_0), 0, cos(pi/2-beta_0)] ;

%vector from CoM to wing hinge
% r_hinge = 0.0006^2 + 0.00018^2 ; %need to check this (in meters). distance from CoM to hinge
% % hinge_vec = R_sp*(1e-3)*[0.18 ; 0; 0.6] ; %body coordinates
% hinge_comp_vert = (1e-3)*0.6 ;
% hinge_comp_fwd = (1e-3)*0.18 ; 
% hinge_vec = R_sp*[hinge_comp_fwd ; 0; hinge_comp_vert] ; %body coordinates

% r_hinge = 0.0003 ; %(1e-3)*0.35 ;%(1e-3)*0.4 ; % meters % from Ristroph, et al., 2013
hinge_vec = [r_hinge; 0; 0] ; % now we're assuming that the vector from body CM to hinge is purely along long body axis
% -------------------------------------------
%% parameters for idealized wing kinematics
psi_m = 57.5 ; % 62.7 ; %60 ;  %deg. This the amplitude for wing pitch %53
psi_0 = 90.0 ; % 90.6 ;  %deg. midpoint of wing pitch % 90
del_psi = -85.0*(pi/180) ; % -76.2*(pi/180) ; %72.4*(pi/180) ; %rad. Wang uses -72.4 degrees (wing pitch phase offset)
C = 2.4 ; % shape factor for wing pitch (aka rotation) angle. C > 0

phi_0 = 97.3 ; %97.5 ; %97.5 ;  %deg. midstroke angle %90 % 100 % 97.7580 ; %
phi_m = 67.6 ; % 67.5 ; %68.2 ;  %deg %63 %70 % 63.5201 ; %
phi_b = phi_0 + phi_m ;  % dorsal limit of stroke angle (units as phi_0, phi_m)
phi_f = phi_0 - phi_m ;  % ventral limit of stroke angle (units as phi_0, phi_m)
K = .7 ;  % shape factor for stroke angle. 0 <= K <= 1

theta_0 = 0 ;  %deg, midpoint of deviation angle
theta_m = 0 ;   %deg, amplitude of deviation angle
del_theta = 0 ;  %rad, phase offset of deviation angle

% ---------------------------------------------
%% lift and drag coefficient params
% axis of wing pitching rotation + rotational force coeff (see Sane & Dickinson, 2002)
x_hat_0 = 0.25 ; 
C_rot_theo = pi*(0.75 - x_hat_0) ; 

% lift and drag coefficients (for approximate form, see Whitney & Wood,
% 2010)
C_L_max = 1.8 ; 
C_D_0 = 0.4 ; 
C_D_max = 3.4 ; 

% ----------------------------------------------------
%% BODY drag coefficients (from Ristroph et al, 2013)
% T_F = 1/drag_beta = 14*5.1 (ms)
drag_beta = 14.0 ; % Hz 
% g * drag_gamma / drag_beta = M*g*h / MOI (where h is vertical distance from CM to wing hinge)
drag_gamma = drag_beta*((body_mass * hinge_vec(3)) / MOI_pitch); % 1/(m*s)
% drag_delta * drag_gamma / drag_beta = M*h^2*drag_beta / MOI
drag_delta = (body_mass * hinge_vec(3)^2 * drag_beta^2)/(MOI_pitch * drag_gamma) ; % m/s
% T_p = 1/drag_epsilon = 14*6.1 (ms)
drag_epsilon = 11.7 ; % Hz 

% ----------------------------------------------------
%% BODY drag coefficients (from Cheng et al, 2009)
C_friction = 0.52e-12 ; % N*m*s

% ----------------------------------------------------
%% perturbation characteristics
pulseStart = 0.0 ; %seconds
pulseEnd = 0.007 ; % seconds
pulseStrength = 7e3 ; %7e3 ; %rad/s^2 (so it's an acceleration not a torque). Roughly estimated from movies %7e3

% -----------------------------------------------------
%% default simulation PI controller values
K_i_sim = 0.437556 ;        % integral gain (unitless)
K_p_sim = 0.005728 ;        % prop gain (sec)
deltaT_sim = 0.004296 ;     %time delay (sec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add values to structure
qs_params = struct() ; 

% -----------------------------------
% physical constants
qs_params.g = g ; 
qs_params.airDensity = airDensity ; 

% -----------------------------------
% computational parameters
qs_params.fitType = fitType ; 
qs_params.bodyVelSmoothWin  = bodyVelSmoothWin ; 
qs_params.wingKinDiffType = wingKinDiffType ; 

% -----------------------------------
% wing morphological parameters
qs_params.freq = freq ; 
qs_params.omega = omega ;
qs_params.span = span ; 
qs_params.chord = chord ; 
qs_params.S = S ; 
qs_params.wing_area = S ; 

qs_params.r22_S = r22_S ; % note: r_XY means subscript X, superscript Y
qs_params.r12_S = r12_S ; 
qs_params.r1_S = r1_S ; 
qs_params.r2_S = r2_S ; 

qs_params.p_beta  = p_beta ; 
qs_params.q_beta  = q_beta ; 
qs_params.Beta_pq = Beta_pq ;

qs_params.c2_r_int = c2_r_int ; % this is integral [0,1] of c_hat^2*r_hat d(r_hat)
qs_params.c2_int = c2_int ; % this is integral [0,1] of c_hat^2 d(r_hat)

qs_params.mean_chord = mean_chord ; 
% -----------------------------------
% body morphological parameters
qs_params.body_mass = body_mass ; 
qs_params.body_length = body_length ; 
qs_params.thorax_width = thorax_width ; 

qs_params.MOI_pitch     = MOI_pitch ; 
qs_params.MOI_roll      = MOI_roll ; 
qs_params.MOI_yaw       = MOI_yaw ; 
qs_params.MOI_rollYaw   = MOI_rollYaw ; 

qs_params.Ixx           = Ixx ; 
qs_params.Iyy           = Iyy ; 
qs_params.Izz           = Izz ; 
qs_params.Ixz           = Ixz ; 

qs_params.Ixx_sp        = Ixx_strokePlane ; 
qs_params.Iyy_sp        = Iyy_strokePlane ; 
qs_params.Izz_sp        = Izz_strokePlane ; 

qs_params.beta_0 = beta_0 ; 

% -----------------------------------
% stroke plane and wing base
qs_params.sp_angle = beta_0 ; 
qs_params.R_sp = R_sp ;
qs_params.r_hinge = r_hinge ; 
qs_params.hinge_vec = hinge_vec ; 

% -----------------------------------
% idealized kinematics
qs_params.psi_m = psi_m ; 
qs_params.psi_0 = psi_0 ; 
qs_params.del_psi = del_psi ;
qs_params.C = C ; 

qs_params.phi_m = phi_m ; 
qs_params.phi_0 = phi_0 ; 
qs_params.phi_f = phi_f ; 
qs_params.phi_b = phi_b ; 
qs_params.K = K ; 

qs_params.theta_m = theta_m ; 
qs_params.theta_0 = theta_0 ; 
qs_params.del_theta = del_theta ;

% ------------------------------------
% force coefficients
qs_params.x_hat_0 = x_hat_0 ; 
qs_params.C_rot_theo = C_rot_theo ; 

qs_params.C_L_max = C_L_max ; 
qs_params.C_D_0 = C_D_0 ; 
qs_params.C_D_max = C_D_max ; 

% ------------------------------------
% BODY drag coefficients
% from Leif:
qs_params.drag_beta     = drag_beta ; 
qs_params.drag_gamma    = drag_gamma ; 
qs_params.drag_delta    = drag_delta ; 
qs_params.drag_epsilon  = drag_epsilon ; 

% from Cheng:
qs_params.C_friction = C_friction ;

% ------------------------------------
% perturbation characteristics
qs_params.pulseStart    = pulseStart ; 
qs_params.pulseEnd      = pulseEnd ; 
qs_params.pulseStrength = pulseStrength ; 

% ------------------------------------
% default simultion controller params
qs_params.K_i_sim       = K_i_sim ; 
qs_params.K_p_sim       = K_p_sim ; 
qs_params.deltaT_sim    = deltaT_sim ;

end

