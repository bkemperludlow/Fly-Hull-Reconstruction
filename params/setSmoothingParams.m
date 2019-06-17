%--------------------------------------------------------------------------
% function to generate parameter structure that will store all of the
% information about angle and position smoothing in fly analysis
%--------------------------------------------------------------------------
function params = setSmoothingParams()

%% body cm 
body_cm_filt_lvl = 100 ; % Hz, for low-pass filter
body_cm_filt_order = 3 ; 
body_cm_est_err = 0.25 ; % in voxels (for spline)

body_sGolayOrder = 4 ;    % poly order for savitzky golay
body_smoothWindow = 401 ; % smoothing window for savitzky golay
spanWindow = 0.2 ; % regression span for loess

%% body yaw
yaw_filt_lvl = 100 ; % Hz, for low-pass filter

%% body roll
roll_filt_lvl = 200 ; % Hz, for low-pass filter
roll_est_err = 0.025 ; % 0.25 ; % deg, for spline fit

%% body pitch
pitch_filt_lvl = 100 ; % Hz, for low-pass filter

%% wing stroke
phi_est_err = 1.75 ; % degrees, estimated error for spline fit
phi_filt_lvl = 500 ; % Hz, for low-pass filter

phi_sgolay_k = 7 ; % polynomial order for savitzky golay filt
phi_sgolay_f = 21 ; % window size for savitzky golay filt
%% wing elevation
theta_est_err = 1.25 ; % degrees, estimated error for spline fit
theta_filt_lvl = 1000 ; % Hz, for low-pass filter

theta_sgolay_k = 5 ; % polynomial order for savitzky golay filt
theta_sgolay_f = 11 ; % window size for savitzky golay filt

%% wing pitch
eta_est_err = 2 ; % degrees, estimated error for spline fit
eta_filt_lvl = 1000 ; % Hz, for low-pass filter

eta_sgolay_k = 5 ; % polynomial order for savitzky golay filt
eta_sgolay_f = 11 ; % window size for savitzky golay filt

%% general wing
wing_hampel_k = 7 ; % hampel filter window size 
wing_hampel_nsigma = 2.5 ; % numbr of std deviations for hampel filter
wing_point_filt_lvl = 2000 ; % Hz, low pass filter level
%--------------------------------------------------------------------------
%% put values in params struct
%--------------------------------------------------------------------------
params.body_cm_filt_lvl = body_cm_filt_lvl ; 
params.body_cm_filt_order = body_cm_filt_order ;
params.body_cm_est_err = body_cm_est_err ;

params.body_sGolayOrder = body_sGolayOrder ;    
params.body_smoothWindow = body_smoothWindow; 
params.body_spanWindow = spanWindow ;
%---------------------------------------------
params.yaw_filt_lvl = yaw_filt_lvl ;

%---------------------------------------------
params.roll_filt_lvl = roll_filt_lvl ;
params.roll_est_err = roll_est_err ; 

%---------------------------------------------
params.pitch_filt_lvl = pitch_filt_lvl ;

%---------------------------------------------
params.phi_est_err = phi_est_err ; 
params.phi_filt_lvl = phi_filt_lvl ;

params.phi_sgolay_k = phi_sgolay_k ;
params.phi_sgolay_f = phi_sgolay_f ; 
%---------------------------------------------
params.theta_est_err = theta_est_err ; 
params.theta_filt_lvl = theta_filt_lvl ;

params.theta_sgolay_k = theta_sgolay_k ; 
params.theta_sgolay_f = theta_sgolay_f ; 
%---------------------------------------------
params.eta_est_err = eta_est_err ; 
params.eta_filt_lvl = eta_filt_lvl ;

params.eta_sgolay_k = eta_sgolay_k ; 
params.eta_sgolay_f = eta_sgolay_f ; 
%--------------------------------------------
params.wing_hampel_k = wing_hampel_k ; 
params.wing_hampel_nsigma = wing_hampel_nsigma ; 
params.wing_point_filt_lvl = wing_point_filt_lvl ; 


end