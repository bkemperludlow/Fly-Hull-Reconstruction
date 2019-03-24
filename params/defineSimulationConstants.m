function [freq, omega, span, chord, body_mass, body_length, MOI, g, airDensity,...
    S, r22_S, sp_angle, R_sp, r_hinge, hinge_vec, psi_m, psi_0, delta_psi, C,...
    phi_0, phi_m, K, theta_0, theta_m, delta_theta, beta_0 ] = defineSimulationConstants()

freq = 250 ; %Hz %250
omega = 2*pi*freq ; 
span = .0025 ; %meters %.002
chord = .0007 ; %meters %.001
body_mass = 1.1e-6 ; %kg
body_length = .0024 ; %meters
%body_width = .0012 ; %meters
MOI = .506e-12 ; %N*m*s^2
beta_0 = pi/4 ; %resting pitch angle, radians

g = 9.8 ; %m/s^2 
airDensity = 1.2041 ; %kg/m^3
S = pi*(span/2)*(chord/2) ; %area of ellipse
r22_S = .313 ; %Sane+Dickinson (2002) gives .4, Cheng et al. gives .313

sp_angle = pi/4 ; %RADIANS. angle between body axis and stroke plane (assumed fixed). Using     
%rotates the velocities and such into the stroke plane in the body axis
R_sp = [ 1, 0, 0 ; 0, cos(pi/2-sp_angle), (-1)*sin(pi/2-sp_angle); ...
            0, sin(pi/2-sp_angle), cos(pi/2-sp_angle)] ;

%vector from CoM to wing hinge
r_hinge = .0006^2 + .00018^2 ; %need to check this (in meters). distance from CoM to hinge
hinge_vec = R_sp*[0; .00018; .0006] ; %body coordinates

% Wing kinematics
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

