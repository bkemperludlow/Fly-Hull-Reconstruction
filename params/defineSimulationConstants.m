% -------------------------------------------------------------------------
% Generic parameters for calculating aerodynamic forces/torques and running
% flight simulations
%
% -------------------------------------------------------------------------
function [freq, omega, span, chord, body_mass, body_length, MOI, g, ...
    airDensity,S, r22_S, sp_angle, R_sp, r_hinge, hinge_vec, psi_m, ...
    psi_0, delta_psi, C, phi_0, phi_m, K, theta_0, theta_m, delta_theta,...
    beta_0, r1_S, c2_r_int, c2_int, C_rot_theo] = defineSimulationConstants()
% ------------------------------------
%% physical constants
g = 9.8 ; %m/s^2 
airDensity = 1.2041 ; %kg/m^3

% ------------------------------------
%% body morphological parameteres
body_mass = 1.1e-6 ; %kg
body_length = .0024 ; %meters
%body_width = .0012 ; %meters
MOI = .506e-12 ; %N*m*s^2
beta_0 = pi/4 ; %resting pitch angle, radians

% ------------------------------------
%% wing morphological parameters
freq = 250 ; %Hz %250
omega = 2*pi*freq ; 
span = .0025 ; %meters %.002
chord = .0007 ; %meters %.001
S = pi*(span/2)*(chord/2) ; %area of ellipse
r22_S = .313 ; %2nd moment of wing inertia. Sane+Dickinson (2002) gives .4, Cheng et al. (2009) gives .313

% can derive other wing shape parameters based on r22_S (see Ellington
% (1984) and  Whitney & Wood (2010))
r2_S = sqrt(r22_S) ; % non-dimensional radius of 2nd moment of wing area
r1_S = (r2_S/0.929).^(1/0.732) ; % non-dimensional radius of 1st moment of wing area
r12_S = r1_S.^2 ; % squared 1st moment of wing area
% parameters for beta distribution describing wing shape
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
% ------------------------------------
%% stroke plane and wing base
sp_angle = pi/4 ; %RADIANS. angle between body axis and stroke plane (assumed fixed). Using     
%rotates the velocities and such into the stroke plane in the body axis
R_sp = [ 1, 0, 0 ; 0, cos(pi/2-sp_angle), (-1)*sin(pi/2-sp_angle); ...
            0, sin(pi/2-sp_angle), cos(pi/2-sp_angle)] ;

%vector from CoM to wing hinge
r_hinge = .0006^2 + .00018^2 ; %need to check this (in meters). distance from CoM to hinge
hinge_vec = R_sp*[0; .00018; .0006] ; %body coordinates

% axis of wing pitching rotation (see Sane & Dickinson, 2002)
x_hat_0 = 0.25 ; 
C_rot_theo = pi*(0.75 - x_hat_0) ; 
% -------------------------------------------
%% parameters for idealized wing kinematics
psi_m = 53 ; %deg. This the amplitude for wing pitch
psi_0 = 90 ; 
delta_psi = 72.4*(pi/180) ; %rad. Wang uses -72.4 degrees
C = 2.4 ;

phi_0 = 95; %deg. midstroke angle %90
phi_m = 75 ; %deg %63
K = .7 ;

theta_0 = 0; %deg
theta_m = 0 ; 
delta_theta = 0 ; 

%Control parameters (Dt defined in main script)
%K_i = .3 ;
%K_p = .007 ; %.01 ;

%Perturbation characteristics
%pulseStart = .01 ; %seconds
%pulseEnd = .015 ;
%pulseStrength = 7e3 ; %rad/s^2 (so it's an acceleration not a torque). Roughly estimated from movies


end

