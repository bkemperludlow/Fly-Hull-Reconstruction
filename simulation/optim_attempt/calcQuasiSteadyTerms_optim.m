% -------------------------------------------------------------------------
% calculate quasi-steady force terms (translational, rotational, and added
% mass)
%
% INPUTS:
%
%
% OUTPUTS:
%   - F_t: translational force
%   - F_rot: rotational force
%   - F_a: added-mass force
% -------------------------------------------------------------------------
function [F_tot, F_t, F_rot, Lhat, Dhat, wingNormalHat] = ...
    calcQuasiSteadyTerms_optim(U_t, spanHat, chordHat, omegaWingWing, ...
        alpha)
% ----------------------------
%% inputs and params
% addedMassFlag = false ; % calculate added mass force?
% smoothFlag = false ; % smooth wing tip velocity?

% read out params from struct
airDensity = 1.2041 ; %kg/m^3
span = 0.002125 ; %in meters (was .002)
mean_chord = 0.0007/1.2681 ; % in meters
S = 5.7375e-06 ; % wing area in m^2
r22_S = 0.313 ; % second moment of wing area
beta_0 = pi/4 ; % resting pitch angle in radians

C_rot_theo = 0.25*pi  ; % coefficient of rotational force (Sane & Dickinson, 2002)
c2_r_int = 0.5395 ; % area moment of some kind 
% c2_int = params.c2_int ; % area moment of some kind

% matrix to pitch down by thetaB0 (~45deg) w.r.t body axis
%M2 = eulerRotationMatrix(0, -1*beta_0, 0) ;
sp2body = eulerRotationMatrix(0, beta_0, 0) ; % stroke plane to body rotation matrix
% ------------------------------------------------
%% translational force
% lift and drag coeffs
C_L_max = 1.8 ; % take hard-coded values
C_D_0 = 0.4 ;
C_D_max = 3.4 ;

% get current lift and drag coefficients for current angle of attack
C_L = C_L_max*sin(2*alpha) ; % lift coeff
C_D = ((C_D_max + C_D_0)/2)*(1 - cos(2*alpha)) ; % drag coeff

% squared magnitude of wing tip velocity
Usq = sum(U_t.*U_t,2) ; 

% % smooth Usq?
% if smoothFlag
%    Usq = smooth(Usq, 3) ;  
% end

% ---------------------
% force magnitudes
fLift = 0.5 * airDensity * S * r22_S * Usq .* C_L ;
fDrag = 0.5 * airDensity * S * r22_S * Usq .* C_D ;

% direction vectors
Uhat = U_t./myNorm(U_t) ; % unit vector of wing tip velocity
Lhat = cross(Uhat, spanHat) ; % perpendicular to Uhat
Lhat = Lhat ./ myNorm(Lhat) ; % but still need to determine sign
% in body fixed frame, should always have positive z component
Lsign = sign(Lhat(:,3)) ;
Lhat((Lsign < 0),:) = -1*Lhat((Lsign < 0),:) ;

Dhat = -Uhat ;

% combine to get full force vectors
F_L = fLift.*Lhat ;
F_D = fDrag.*Dhat ;

F_t = F_L + F_D ;

% ------------------------------------------------
%% rotational force
% wing rotation about span from (Sane & Dickinson, 2002)
omega_rot = omegaWingWing(:,1) ;

% magnitude of rotational force
fRot = C_rot_theo * airDensity * sqrt(Usq) .* abs(omega_rot) .* ...
    mean_chord^2 .* span .* c2_r_int ;

% should act perpendicular to wing surface
wingNormalHat = cross(spanHat, chordHat) ; % need to check sign
wingNormalHat = wingNormalHat./myNorm(wingNormalHat) ;

% one way of ensuring we have the right sign for wing normal is to compare
% it to the total translational force vector--this should also be normal to
% wing, in general
F_t_hat = F_t./myNorm(F_t) ;
wingNormalHat = sign(dot(wingNormalHat, F_t_hat,2)).*wingNormalHat ;

% get full rotational force vector
F_rot = fRot .* wingNormalHat ;

% ------------------------------------------------
% %% added mass force
% if addedMassFlag
%     % added mass force has 2 terms. From (Sane & Dickinson, 2002). Note 
%     % that the "phi" in their notation refers to "angular position" of 
%     % wing, not just stroke angle (I think). So first thing is to get 
%     % (d/dt)omegaBodyWing :
%     c_omega_wb_x = fit(t, omegaWingBody(:,1), fitType) ;
%     c_omega_wb_y = fit(t, omegaWingBody(:,2), fitType) ;
%     c_omega_wb_z = fit(t, omegaWingBody(:,3), fitType) ;
%     
%     omega_wb_x_dot = differentiate(c_omega_wb_x,t) ;
%     omega_wb_y_dot = differentiate(c_omega_wb_y,t) ;
%     omega_wb_z_dot = differentiate(c_omega_wb_z,t) ;
%     omega_wb_dot = [omega_wb_x_dot,omega_wb_y_dot,omega_wb_z_dot] ;
%     
%     % get norm
%     omega_wb_norm = myNorm(omegaWingBody) ;
%     omega_wb_dot_norm = myNorm(omega_wb_dot) ;
%     
%     % also need derivatives of angle of attack w.r.t. time
%     c_alpha = fit(t, alpha, fitType) ;
%     [alpha_dot, alpha_ddot] = differentiate(c_alpha, t) ;
%     
%     % now calculate the two terms in the added mass force (magnitudes)
%     fAddedMass1 = (airDensity * (pi/4) * span^2 * mean_chord^2) .*...
%         c2_r_int .*(omega_wb_dot_norm .* sin(alpha) + ...
%         omega_wb_norm.*alpha_dot.*cos(alpha))  ;
%     fAddedMass2 = -1 * alpha_ddot .* airDensity .* (pi/16) .*...
%         mean_chord^3 .* span .* c2_int ;
%     
%     % sum to get total added-mass force
%     fAddedMass = fAddedMass1 + fAddedMass2 ;
%     
%     % added-mass force also acts normal to wing
%     F_a = fAddedMass .* wingNormalHat ;
% else
%     F_a = zeros(size(F_t)) ;
% end

% ----------------------------------------
%% sum components to get total force
F_tot = F_t + F_rot ; 

% --------------------------------------------------------------------
% %% rotate force terms from stroke plane frame to body frame
F_tot = (sp2body*F_tot')' ; 
F_t = (sp2body*F_t')' ; 
F_rot = (sp2body*F_rot')' ; 
% F_a = (sp2body*F_a')' ; 

% also rotate lift and drag unit vectors
Lhat = (sp2body*Lhat')' ; 
Dhat = (sp2body*Dhat')' ; 

end