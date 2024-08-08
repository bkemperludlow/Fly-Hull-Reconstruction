% -------------------------------------------------------------------------
% function to use of testing the force/torque produced by artificial wing
% kinematics as we vary one parameter
% -------------------------------------------------------------------------
% ------------------------
%% get/set param values
% which parameter will we vary?
varyParam1 = 'delta_phi_f' ; % 'delta_phi_f' ; %'omega' ;
varyParam2 = 'pitch_ang' ; % 'phi_b' % NB: can currently only vary body kin in param2
N_vals = 51 ; % 50

% define params structure
params = defineQuasiSteadyParams() ;

% get wing_kin_params vector, but update with values from x
wing_kin_params = read_wing_kin_params(params) ;

% get body_params vector ;
body_params = read_body_params(params) ;

% other params
thetaB0 = params.beta_0 ; % stroke plane to body coords angle
body_mass = params.body_mass ;
g = params.g ;
thorax_width = params.thorax_width ;

% misc settings
addedMassFlag = false ;
plotFlag = false ;

% ------------------------------------------------
%% define ranges over which we'll test params vals
omega_range = 2*pi*linspace(100, 300, N_vals) ;
phi_f_range = (pi/180)*linspace(-10, 70, N_vals) ;
phi_b_range = (pi/180)*linspace(100, 200, N_vals) ;
delta_phi_f_range = (pi/180)*linspace(-25, 25, N_vals) ;

span_range = (1e-3)*linspace(1, 3, N_vals) ;
hinge_x_range = (1e-3)*linspace(0.01, 0.5, N_vals) ;

% range to use for varying pitch velocity (in rad/s)
% NB: this is a big range
pitchVelMin = -150 ; %-60 ; % rad/s 
pitchVelMax = 150; %  60 ;  % rad/s
pitch_vel_range = linspace(pitchVelMin, pitchVelMax, N_vals) ; 

% range to use for varying pitch angle (in radians)
pitchAngMin = (pi/180)*0 ;  % radians
pitchAngMax = (pi/180)*90 ; % radians
pitch_ang_range = linspace(pitchAngMin, pitchAngMax, N_vals) ; 
% ------------------------------------------------
%% initialize output matrices
F_x_all = nan(N_vals, N_vals) ;
F_z_all = nan(N_vals, N_vals) ;
T_y_all = nan(N_vals, N_vals) ;

% ------------------------------------------------
%% loop over 1st param values
for k = 1:N_vals
    % get current value of given param
    switch varyParam1
        case 'omega'
            wing_kin_params(1) = omega_range(k) ;
        case 'phi_f'
            wing_kin_params(2) = phi_f_range(k) ;
        case 'phi_b'
            wing_kin_params(3) = phi_b_range(k) ;
        case 'span'
            params = defineQuasiSteadyParams('span',span_range(k)) ;
            body_mass = params.body_mass ;
            body_params(1) = params.span ;
        case 'hinge_x'
            params = defineQuasiSteadyParams('r_hinge',hinge_x_range(k)) ;
            body_params(2) = params.r_hinge ;
        case 'delta_phi_f'
            wing_kin_params(12) = delta_phi_f_range(k) ; 
        otherwise
            % do nothing -- as in the case when we want to vary body pitch
            % velocity
            
    end
    
    % -------------------------------------
    %% loop over 2nd param vals
    for m = 1:N_vals
        % since we might vary pitch angle/velocity, create an array of 
        % zeros for that now (which we can overwrite if we're in the 
        % correct case)
        pitch_vel_curr = 0 ;
        pitch_ang_curr = thetaB0 ;
        % --------------------------------------
        %% get current value of 2nd param
        switch varyParam2
            case 'omega'
                wing_kin_params(1) = omega_range(m) ;
            case 'phi_f'
                wing_kin_params(2) = phi_f_range(m) ;
            case 'phi_b'
                wing_kin_params(3) = phi_b_range(m) ;
            case 'span'
                params = defineQuasiSteadyParams('span',span_range(m),...
                    'r_hinge',params.r_hinge) ;
                body_mass = params.body_mass ;
                body_params(1) = params.span ;
            case 'hinge_x'
                params = defineQuasiSteadyParams('span', params.span, ...
                    'r_hinge',hinge_x_range(m)) ;
                body_params(2) = params.r_hinge ;
            case 'delta_phi_f'
                wing_kin_params(12) = delta_phi_f_range(k) ; 
            case 'pitch_vel'
                % assign body pitch velocity 
                % NB: this is sort of crude, since we're not varying the
                % body angle in accordance with the pitch velocity;
                % however, it should still be a good approximation?
                pitch_vel_curr = pitch_vel_range(m) ; 
            case 'pitch_ang'
                % assign body pitch angle
                pitch_ang_curr = pitch_ang_range(m) ; 
        end
        
        % -------------------------------------------
        %% get time range (one wingbeat)
        omega = wing_kin_params(1) ;
        ti = 0 ;
        tf = (2*pi)/omega ;
        t = linspace(ti, tf, 100) ;
        dt = mean(diff(t)) ; 
        
        % shift so that we start at beginning of downstroke
        t_shift = (pi/2)/omega ;
        t = t + t_shift ;
        
         % -------------------------------------------------------------
        %% body movement -- set to zero unless we're varying pitch vel
        velBodyBody = zeros(length(t),3) ;
        bodyYPR = zeros(length(t),3) ;
        bodyYPR(:,2) = pitch_ang_curr.*ones(length(t),1) ; % set steady pitch (should default to thetaB0 unless we're intentionally adjusting it)
        bodyYPR_dot = zeros(length(t),3) ;
        
        % set pitch velocity (probably zero, unless we're varying it)
        bodyYPR_dot(:,2) = repmat(pitch_vel_curr, length(t),1) ; 
        
%         % try to adjust the body pitch angle accordingly?
%         bodyYPR(:,2) = bodyYPR(1,2) + dt.*cumsum(bodyYPR_dot(:,2)) ; 
        
        % ------------------------------------------------
        %% calculate force and torque
        [F_tot, T_tot] = fakeWingKinForceAndTorque(wing_kin_params, ...
            body_params, t, velBodyBody, bodyYPR, bodyYPR_dot, params,...
            addedMassFlag, plotFlag) ;
        
        % plot results?
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
            
        end
        
        % -----------------------------------------------------------------
        %% normalize force and torque, take mean, and store in array
        F_mean = nanmean(F_tot)./(body_mass*g) ;
        % normalizing by thorax width here is weird, but w/e
        T_mean = nanmean(T_tot)./(body_mass*g*thorax_width) ;
        
%         % also subtract gravity from F
%         thetaB_mean = nanmean(bodyYPR(:,2)) ;
%         rotM = eulerRotationMatrix(0, thetaB_mean, 0) ;
%         g_norm = [0; 0; 1] ;
%         g_norm_rot = (rotM'*g_norm)' ;
%         
%         F_mean = F_mean - g_norm_rot ;
        % rotate force into lab coordinates 
        M_thetaB0 = eulerRotationMatrix(0, nanmean(bodyYPR(:,2)), 0) ;
        F_mean = (M_thetaB0'*F_mean')' ;
        
        % get individual components we care about
        F_x = F_mean(1) ;
        F_z = F_mean(3) - 1 ;
        T_y = T_mean(2) ;
        
        % store in matrix
        F_x_all(k,m) = F_x ;
        F_z_all(k,m) = F_z ;
        T_y_all(k,m) = T_y ;
        
    end
end

% -------------------------------------------------------------
%% plot results
% first get x and y axis data (varyParam1, as rows, is y axis)
switch varyParam1
    case 'omega'
        y = omega_range./(2*pi) ;
        y_label = 'wbf (Hz)' ;
    case 'phi_f'
        y = (180/pi).*phi_f_range ;
        y_label = '\phi_f (deg)' ;
    case 'phi_b'
        y = (180/pi).*phi_b_range ;
        y_label = '\phi_b (deg)' ;
    case 'delta_phi_f'
        y = delta_phi_f_range ;
        y_label = '\Delta\phi_f (rad)' ;
    case 'span'
        y = span_range ;
        y_label = 'span (m)' ;
    case 'hinge_x'
        y = hinge_x_range ;
        y_label = 'hinge x (m)' ;
    otherwise
        keyboard
end
switch varyParam2
    case 'omega'
        x = omega_range./(2*pi) ;
        x_label = 'wbf (Hz)' ;
    case 'phi_f'
        x = (180/pi).*phi_f_range ;
        x_label = '\phi_f (deg)' ;
    case 'phi_b'
        x = (180/pi).*phi_b_range ;
        x_label = '\phi_b (deg)' ;
    case 'delta_phi_f'
        x = delta_phi_f_range ;
        x_label = '\Delta\phi_f (rad)' ;
    case 'span'
        x = span_range ;
        x_label = 'span (m)' ;
    case 'hinge_x'
        x = hinge_x_range ;
        x_label = 'hinge x (m)' ;
    case 'pitch_vel'
        x = pitch_vel_range ; 
        x_label = 'pitch vel. (rad/s)' ; 
    case 'pitch_ang'
        x = pitch_ang_range ; 
        x_label = 'pitch angle (rad)' ; 
end

% ------------------
% F_x
% ------------------
h_x = figure ; 
imagesc([x(1), x(end)], [y(1), y(end)], F_x_all)

% colorbar/map
val_max = nanmax(abs(F_x_all(:))) ; 
caxis(val_max.*[-1, 1])
colormap(brewermap([],'RdBu')) ; 
colorbar

% axis labels
title('F_x')
ylabel(y_label) 
xlabel(x_label)

% ------------------
% F_z
% ------------------
h_z = figure ; 
imagesc([x(1), x(end)], [y(1), y(end)], F_z_all)

% colorbar/map
val_max = nanmax(abs(F_z_all(:))) ; 
caxis(val_max.*[-1, 1])
colormap(brewermap([],'RdBu')) ; 
colorbar

% axis labels
title('F_z')
ylabel(y_label) 
xlabel(x_label)

% ------------------
% T_y
% ------------------
h_y = figure ; 
imagesc([x(1), x(end)], [y(1), y(end)], T_y_all)

% colorbar/map
val_max = nanmax(abs(T_y_all(:))) ; 
caxis(val_max.*[-1, 1])
colormap(brewermap([],'RdBu')) ; 
colorbar

% axis labels
title('T_y')
ylabel(y_label) 
xlabel(x_label)

% ------------------
% T_y SURFACE PLOT
% ------------------
[xx, yy] = meshgrid(x, y) ; 
h_y_surf = figure ; 
hold on

% % make surface plot
% s = surf(xx, yy, T_y_all, 'EdgeColor','none', 'FaceColor','interp') ;
% 
% % set colorbar properties
% val_max = nanmax(abs(T_y_all(:))) ; 
% caxis(val_max.*[-1, 1])
% colormap(brewermap([],'RdBu')) ; 
% colorbar

plot3(xx(:), yy(:), T_y_all(:), 'k.','MarkerFaceColor','none','MarkerSize',8)

% axis labels
zlabel('T_y')
ylabel(y_label) 
xlabel(x_label)

% try to perform fit to surface? (polynomial model)
sf = fit([xx(:), yy(:)], T_y_all(:), 'poly22') ; 
sfit_plot = surf(xx, yy, sf(xx,yy), 'EdgeColor','none', ...
    'FaceColor','interp') ; 

% set colorbar properties
val_max = nanmax(abs(T_y_all(:))) ; 
caxis(val_max.*[-1, 1])
colormap(brewermap([],'RdBu')) ; 
colorbar

% turn on 3d rotation
rotate3d on

% get tangent to plane at x = y = 0 
[fx, fy] = differentiate(sf, [0,0]) ; 

% for param1 = delta_phi_front and param2 = pitch_vel, get
% fx = 0.0018 ; 
% fy = (180/pi)*0.0287 ;  % NB: this fy is for deltaPhiFront in degrees 

% % --------------
% % summed force
% % --------------
% h_sum = figure ; 
% F_sum = F_x_all.^2 + F_z_all.^2 ; 
% imagesc([x(1), x(end)], [y(1), y(end)], F_sum)
% 
% % colorbar/map
% val_max = nanmax(abs(F_sum(:))) ; 
% % caxis(val_max.*[-1, 1])
% % colormap(brewermap([],'RdBu')) ; 
% colorbar
% 
% % axis labels
% title('F_{sum}')
% ylabel(y_label) 
% xlabel(x_label)
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