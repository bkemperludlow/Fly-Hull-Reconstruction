% -------------------------------------------------------------------------
% function to make a movie of 3D flie based on simulation results
%
% INPUTS
%   - sol: structure output of diff eq solver containing solution to flight
%       simulation. See "simulatePitch.m"
%   - params: structure of parameters used for quasi-steady calculations.
%       See "defineQuasiSteadyParams.m"
%
% OUTPUTS:
%   -
% -------------------------------------------------------------------------
function [] = makeSimulationMovie(sol, params, controllerTerms, thetaB0,...
    savePath, prefixStr, controlFlag, pertFlag)
% ---------------------
%% inputs and params
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams() ;
end
if ~exist('controllerTerms','var') || isempty(controllerTerms)
    % individual controller terms
    K_i = 0.437556 ;
    K_p = 0.005728 ;
    deltaT = 0.002148 ;
    
    % add to vector
    controllerTerms = [K_i, K_p, deltaT] ;
end
if ~exist('thetaB0','var') || isempty(thetaB0)
    thetaB0 = params.beta_0 ;
end
if ~exist('savePath','var') || isempty(savePath)
    [savePath, ~, ~] = fileparts(mfilename('fullpath')) ;
    %t_now = now ;
    %     savePath = fullfile(savePath, ['temp_data_' num2str(t_now)]) ;
    savePath = fullfile(savePath, 'temp_data') ;
end
if ~exist('prefixStr','var') || isempty(prefixStr)
    prefixStr = '' ;
end
if ~exist('controlFlag','var') || isempty(controlFlag)
    controlFlag = false ;
end
if ~exist('pertFlag','var') || isempty(pertFlag)
    pertFlag = false ;
end

% make sure savePath exists
if ~exist(savePath,'dir')
    mkdir(savePath)
end

% -------------------------
% parameters for video
vid_profile = 'MPEG-4' ; 
vid_framerate = 80 ; 

switch vid_profile
    case 'MPEG-4'
        vid_file_ext = '.mp4' ; 
    case {'Uncompressed AVI', 'Motion JPEG AVI'}
        vid_file_ext = '.avi' ; 
    otherwise
        fprintf('Invalid video profile: %s \n', vid_profile)
        keyboard
end

% fly size/resolution
scale = 35 ;
flyResolution =  100 ;

% gridlines on fly?
gridFlag = false ; 

% view angles
az = 0; %61 ; %54
el = 0 ; %23 ; % 12

% background color for fly
bg_color = 0.0*[1, 1, 1] ; 

% lighting stuff
fly_material = 'shiny' ;
fly_ambient_light = 0.9 ;
fly_color_scheme = 'sim' ;

% pin type
if pertFlag
    pinType = 2 ; % pitch
else
    pinType = 0 ; % no pin
end
% ------------------------------------------------
%% read out data from diff eq solution struct
% time range of simulation
dt = 0.125e-3 ; % seconds. same time resolution as real data
ti = min(sol.x) ;
tf = max(sol.x) ;
t = ti : dt : tf ;

N_frames = length(t) ;

% get body position and body pitch angle
sint = deval(sol, t) ;
x_body = sint(1,:)' ; % x position (fwd/back) in body frame coordinates
x_dot_body = sint(2,:)' ; % x vel (fwd/back) in body frame coordinates
z_body = sint(3,:)' ; % z position (up/down) in body frame coords
z_dot_body = sint(4,:)' ; % z vel (fwd/back) in body frame coordinates
thetaB = sint(5,:)' ; % body pitch, radians
thetaB_dot = sint(6,:)' ; % body pitch velocity , rad/s

% ---------------------------------------------------------------------
% convert body position from body coords to lab coords by integrating
% velocity
x_dot_lab = zeros(N_frames,1) ; 
z_dot_lab = zeros(N_frames,1) ; 

% loop over frames and rotate
for k = 1:N_frames
    rotM = eulerRotationMatrix(0, thetaB(k),0) ;
    
    % rotate velocity vector
    vel_vec = [x_dot_body(k) ; 0 ; z_dot_body(k)] ;
    vel_vec_lab = rotM'*vel_vec ;
    x_dot_lab(k) = vel_vec_lab(1) ;
    z_dot_lab(k) = vel_vec_lab(3) ;
    
end

% perform integration
x_lab = dt.*cumtrapz(x_dot_lab) ;
z_lab = dt.*cumtrapz(z_dot_lab) ; 

% ----------------------------------------------------------------
% compile into matrices
bodyCM = [x_lab, zeros(N_frames,1), z_lab] ;
bodyYPR = [zeros(N_frames,1), thetaB, zeros(N_frames,1)] ;

% scale body position by overall scale of fly drawing
% extra factor of 1000 for meter -> millimeter conversion
% extra factor of 2 because thorax radius is about 2x too big in drawFly3D
bodyCM = (1000*scale).*bodyCM ;  %(2*1000*scale).*bodyCM ; 

% get axis limits from body CM locations
pad_size = 5*scale ; 
xlim = [min(bodyCM(:,1)) - 3*pad_size, max(bodyCM(:,1)) + 3*pad_size] ;
ylim = [min(bodyCM(:,2)) - pad_size, max(bodyCM(:,2)) + pad_size] ;
zlim = [min(bodyCM(:,3)) - pad_size, max(bodyCM(:,3)) + pad_size] ;
axlim = [xlim, ylim, zlim] ; 

% ------------------------------------------
%% get wing kinematics over time range
% general wing params
wing_kin_params = read_wing_kin_params(params) ;

% if we were using control in the simulation, need to get delta_phi_f for
% each frame
if controlFlag
    % read out controller coeffs
    K_i = controllerTerms(1) ;
    K_p = controllerTerms(2) ;
    deltaT = controllerTerms(3) ;
    
    % initialize arrays to store angle data
    wingAngleMatR = zeros(N_frames, 3) ;
    % wingAngleMatL = zeros(N_frames, 3) ;
    
    % loop over time values
    for k = 1:N_frames
        % get lagged values (if no measured values exist, use steady vals)
        if (t(k) < deltaT)
            thetaB_lag = thetaB0 ;
            thetaB_dot_lag = 0 ;
        else
            t_lag = t(k) - deltaT ;
            [~, lag_ind] = min(abs(t - t_lag)) ;
            thetaB_lag = thetaB(lag_ind) ;
            thetaB_dot_lag = thetaB_dot(lag_ind) ;
        end
        
        % get output of PI controller
        delta_phi_f = pitch_controller_func(thetaB_lag, thetaB_dot_lag, ...
            thetaB0, K_i, K_p) ;
        
        % add to wing_kin_params and get current wing kin values
        wing_kin_params_curr = wing_kin_params ;
        wing_kin_params_curr(12) = delta_phi_f ;
        
        [wingAngleMatR(k,:), ~] = getFakeWingKin(wing_kin_params_curr, ...
            'R', t(k), false) ;
    end
    
    % don't actually need sign flip on right wing stroke angle, so remove
    % that and make left wing a copy of right
    wingAngleMatR(:,1) = -1.*wingAngleMatR(:,1) ;
    wingAngleMatL = wingAngleMatR ;
else
    % if not using control, this should be simple
    [wingAngleMatR, ~] = getFakeWingKin(wing_kin_params, 'R', t, false) ;
    
    % don't actually need sign flip on right wing stroke angle, so remove
    % that and make left wing a copy of right
    wingAngleMatR(:,1) = -1.*wingAngleMatR(:,1) ;
    wingAngleMatL = wingAngleMatR ;
end

% ------------------------------------------------------------------
%% initialize video writer object
% % video filename
vid_filename = fullfile(savePath, [prefixStr 'simulation' vid_file_ext]) ;

% video writer object
vid_writer = VideoWriter(vid_filename, vid_profile ) ; 

% set video params
vid_writer.FrameRate = vid_framerate ; 

% ------------------------------------------------------------------
%% loop over frames and draw fly in current position
% initialize figure and axis
h_main = figure('units','normalized','outerposition',[0 0 1 1]) ; % full screen
ax = axes(h_main) ;

% draw 3D fly
[flyGrp, ~, rightWingGrp, leftWingGrp, dL, ~, ~] = draw3Dfly(ax, ...
    scale, flyResolution, pinType, thetaB0, gridFlag, fly_color_scheme );

% fly material properties
material(flyGrp, fly_material)
c = findobj(flyGrp,'Type','surface');
set(c,'ambientstrength',fly_ambient_light); % note that this line comes after "material shiny"

% set view 
view(az,el)

% set axis limits, make axis equal, and remove lines/ticks
axis(ax, axlim, 'equal', 'off','manual')

% set color background
set(ax, 'Color', bg_color)
set(h_main, 'Color', bg_color)

% --------------------------
% open video writer
% --------------------------
open(vid_writer)

% --------------------------
% loop over frames
for ind = 1:N_frames
    % set fly DOF
    setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, bodyCM(ind,:), ...
        bodyYPR(ind,:), wingAngleMatR(ind,:), wingAngleMatL(ind,:),...
        thetaB0, dL) ;
    
    axis(ax, axlim, 'equal', 'off','manual')
    
    % grab frame 
    frame = getframe(ax) ; 
    
%     % pad with zeros if frame dimensions aren't even
%     if mod(size(frame.cdata,1),2) ~= 0
%         frame.cdata = padarray(frame.cdata,1,'post') ;
%     end
%     if mod(size(frame.cdata,2),2) ~= 0
%         frame.cdata = padarray(frame.cdata,2,'post') ;
%     end
    % add frame to video writer
    writeVideo(vid_writer, frame) ; 
    
    % print update
    fprintf('Completed frame %d / %d \n', ind, N_frames)
end

% --------------------------
% close video writer
% --------------------------
close(vid_writer)
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

% --------------------------------------------------------------------
%% PI controller function
function delta_phi_f = ...
    pitch_controller_func(thetaB_lag, thetaB_dot_lag, thetaB0, K_i, K_p)

% change in fwd stroke amplitude is determined by body pitch angle and vel
delta_phi_f = K_i.*(thetaB_lag - thetaB0) + K_p.*(thetaB_dot_lag) ;
delta_phi_f = nanmean(delta_phi_f) ;

if abs(delta_phi_f) > 0.8727
    delta_phi_f = sign(delta_phi_f)*0.8727 ;
end

end
