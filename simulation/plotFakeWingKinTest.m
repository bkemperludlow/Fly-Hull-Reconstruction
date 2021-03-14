% -------------------------------------------------------------------------
% script to plot artificial wing kinematics and force/torque just to check
% that things are working okay
% -------------------------------------------------------------------------
%% define QS parameters
params = defineQuasiSteadyParams ;
addedMassFlag = false ;

wingSideList = {'R','L'} ;
wingSignList = [-1, 1] ;

% get time range (one wingbeat)
omega = params.omega ;
ti = 0 ;
tf = (2*pi)/omega ;
t = linspace(ti, tf, 100) ;

% shift so that we start at beginning of downstroke
t_shift = (pi/2)/omega ;
t = t + t_shift ;

% normalize to wingbeat time
t_norm = linspace(0, 1, length(t)) ; 

% get wing_kin_params vector, but update with values from x
wing_kin_params = read_wing_kin_params(params) ;

% get body_params vector ;
body_params = read_body_params(params) ;

% other params
thetaB0 = params.beta_0 ; % stroke plane to body coords angle
body_mass = params.body_mass ;
g = params.g ;
thorax_width = params.thorax_width ;
span = params.span ;

% ------------------------------------------------
%% no body movement here, so set to zero
velBodyBody = zeros(length(t),3) ;
bodyYPR = zeros(length(t),3) ;
bodyYPR(:,2) = thetaB0.*ones(length(t),1) ; % set steady pitch
bodyYPR_dot = zeros(length(t),3) ;

% -------------------------------------------------------
%% get wing kinematics
wingAngleMat = zeros(length(t), 3, 2) ;
for k = 1:length(wingSideList)
    % current wing side (right vs left)
    wingSide = wingSideList{k} ;
    wingSign = wingSignList(k) ;
    
    % ----------------------------
    % get current wing kin
    [wingAngleMat(:,:,k), ~] = getFakeWingKin(wing_kin_params, ...
        wingSide, t, false) ;
end
% --------------------------------------------------
%% plot wing kin
h_kin = figure ;
kin_labels = {'\phi', '\theta', '\eta'} ;
wing_sign_vec = [-1, 1, 1] ;
for dim = 1:3
    subplot(3,1,dim)
    hold on
    plot(t_norm,...
        wing_sign_vec(dim)*(180/pi).*squeeze(wingAngleMat(:,dim,1)),'r')
    if dim == 3
        plot(t_norm, (180/pi).*(pi - squeeze(wingAngleMat(:,dim,2))),'b')
    else
        plot(t_norm, (180/pi).*squeeze(wingAngleMat(:,dim,2)),'b')
    end
    
    % axis properties
    axis tight
    xlabel('Time (wb)')
    ylabel(sprintf('%s (deg)',kin_labels{dim}))
end
% ------------------------------------------------
%% calculate force and torque
[F_tot, T_tot] = fakeWingKinForceAndTorque(wing_kin_params, ...
    body_params, t, velBodyBody, bodyYPR, bodyYPR_dot, params, ...
    addedMassFlag, false) ;

% --------------------------------------------------
%% plot force/torque output (BODY FRAME)
h_force = figure ;
dim_labels = {'x', 'y', 'z'} ;
for dim = 1:3
    subplot(3,1,dim)
    hold on
    plot(t_norm, F_tot(:,dim)./(body_mass*g))
    plot([t_norm(1), t_norm(end)], ...
        (nanmean(F_tot(:,dim))/(body_mass*g)).*[1, 1],'--')
    
    % axis properties
    axis tight
    xlabel('Time (wb)')
    ylabel(sprintf('F_{%s}/mg',dim_labels{dim}))
end

h_torque = figure ;
dim_labels = {'x', 'y', 'z'} ;
for dim = 1:3
    subplot(3,1,dim)
    hold on
    plot(t_norm, T_tot(:,dim)./(body_mass*g*span))
    plot([t_norm(1), t_norm(end)], ...
        (nanmean(T_tot(:,dim))/(body_mass*g*span)).*[1, 1],'--')
    
    % axis properties
    axis tight
    xlabel('Time (wb)')
    ylabel(sprintf('T_{%s}/mgl',dim_labels{dim}))
end

% ----------------------------------------------------
%% convert forces to lab frame
F_tot_lab = zeros(size(F_tot)) ; 
for ind = 1:length(t)
   lab2body = eulerRotationMatrix(bodyYPR(ind,1), bodyYPR(ind,2), ...
       bodyYPR(ind,3)) ; 
   F_tot_lab(ind,:) = (lab2body')*F_tot(ind,:)' ; 
end

% --------------------------------------------------
%% plot force/torque output (LAB FRAME)
h_force_lab = figure ;
dim_labels = {'x', 'y', 'z'} ;
for dim = 1:3
    subplot(3,1,dim)
    hold on
    plot(t_norm, F_tot_lab(:,dim)./(body_mass*g))
    plot([t_norm(1), t_norm(end)], ...
        (nanmean(F_tot_lab(:,dim))/(body_mass*g)).*[1, 1],'--')
    
    % axis properties
    axis tight
    xlabel('Time (wb)')
    ylabel(sprintf('F_{%s}/mg (lab)',dim_labels{dim}))
end