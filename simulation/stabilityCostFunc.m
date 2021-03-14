% -------------------------------------------------------------------------
% cost function to be used to find kinematic/morphological parameters that
% minimize force/torque in longitudinal flight under hover conditions
%
% here we're solving for x, which is a vector of kinematic/morphological
% params of the form:
%   x = [omega, phi_f, phi_b, span, hinge_x] ; 
% -------------------------------------------------------------------------
function [cost, g, H] = stabilityCostFunc(x)
% ------------------------
%% misc params
addedMassFlag = false ; 
plotFlag = false ; 

% --------------------
%% read values from x
omega = x(1) ; 
phi_f = x(2) ; 
phi_b = x(3) ; 
psi_m = x(4) ; 
psi_0 = x(5) ; 
del_psi = x(6) ;
span = x(7) ; 
hinge_x = x(8) ; 

% define params structure
params = defineQuasiSteadyParams('span',span,'r_hinge',hinge_x) ; 

% get time range (one wingbeat)
ti =  0 ; 
tf = (2*pi)/omega ; 
t = (linspace(ti, tf, 100))' ; 

% shift so that we start at beginning of downstroke
t_shift = (pi/2)/omega ; 
t = t + t_shift ; 

% get wing_kin_params vector, but update with values from x
wing_kin_params = read_wing_kin_params(params) ; 
wing_kin_params(1) = omega ; 
wing_kin_params(2) = phi_f ; 
wing_kin_params(3) = phi_b ; 
wing_kin_params(9) = psi_m ; 
wing_kin_params(8) = psi_0 ; 
wing_kin_params(10) = del_psi ; 
% get body_params vector ; 
body_params = read_body_params(params) ; 

% other params
thetaB0 = params.beta_0 ; % stroke plane to body coords angle
body_mass = params.body_mass ; 
g = params.g ; 
thorax_width = params.thorax_width ; 
% span = params.span ; 
% ------------------------------------------------
%% no body movement here, so set to zero
velBodyBody = zeros(length(t),3) ; 
bodyYPR = zeros(length(t),3) ; 
bodyYPR(:,2) = thetaB0.*ones(length(t),1) ; % set steady pitch 
bodyYPR_dot = zeros(length(t),3) ; 

% ------------------------------------------------
%% calculate force and torque
% [F_tot, T_tot] = fakeWingKinForceAndTorque(wing_kin_params, ...
%     body_params, t, velBodyBody, bodyYPR, bodyYPR_dot, params, ...
%     addedMassFlag, plotFlag) ; 

[F_tot, T_tot] = fakeWingKinForceAndTorque_optim(wing_kin_params, ...
    body_params, t, velBodyBody, bodyYPR_dot, params) ; 

if (0)
   M_thetaB0 = eulerRotationMatrix(0, nanmean(bodyYPR(:,2)), 0) ; 
   F_temp = (M_thetaB0'*F_tot')' ;
   F_temp = F_temp./(body_mass*g) ;
   T_temp = T_tot./(body_mass*g*span) ;
   figure ; 
   subplot(3,1,1)
   hold on
   plot(t, F_temp(:,1))
   plot([t(1), t(end)], nanmean(F_temp(:,1)).*[1,1],'--') 
   set(gca,'xlim',[t(1), t(end)])
   
   subplot(3,1,2)
   hold on
   plot(t, F_temp(:,3))
   plot([t(1), t(end)], nanmean(F_temp(:,3)).*[1,1],'--') 
   set(gca,'xlim',[t(1), t(end)])
   
   subplot(3,1,3)
   hold on
   plot(t, T_temp(:,2))
   plot([t(1), t(end)], nanmean(T_temp(:,2)).*[1,1],'--') 
   set(gca,'xlim',[t(1), t(end)])
   
   keyboard
end
% -------------------------------------------------------------------
%% cost is the squared sum of F_x, F_z, T_y, averaged over wingbeat
F_mean = nanmean(F_tot)./(body_mass*g) ; 
% normalizing by thorax width here is weird, but just need F and T to be on
% similar scales so one doesn't dominate the cost function, and the normal
% choices of length scales are all affected by span, which is something
% we're solving for
T_mean = nanmean(T_tot)./(body_mass*g*thorax_width) ; % (body_mass*g*thorax_width) 

% also subtract gravity from F
thetaB_mean = nanmean(bodyYPR(:,2)) ; 
rotM = eulerRotationMatrix(0, thetaB_mean, 0) ; 
% g_norm = [0; 0; 1] ; 
% g_norm_rot = (rotM'*g_norm)' ; 

F_mean = (rotM'*F_mean')' ;  

% get individual components we care about
F_x = F_mean(1) ; 
F_z = F_mean(3) - 1 ; 
T_y = T_mean(2) ; 

% % normalize
% F_x_norm = F_x/(body_mass*g) ; 
% F_z_norm = F_z/(body_mass*g) ; 
% T_y_norm = T_y/(body_mass*g*span) ; 

% get squared sum
cost = sum([F_x, F_z, T_y].^2) ; 
% cost = T_y^2 ;
% cost = F_x^2 ; 
% ----------------------------------------------------------------
%% don't have a good sense for gradient or Hessian, so return 0
if (nargout > 1)
    g = zeros(size(x)) ;
    if (nargout > 2)
        H = zeros(length(x),length(x)) ;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------
%% read out data from params to a wing_kin_params vector
function wing_kin_params = read_wing_kin_params(params)
% -----------------------------------------
% fix delta phi front in this case, since we're explicitly varying phi_f
delta_phi_f = 0 ; 

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