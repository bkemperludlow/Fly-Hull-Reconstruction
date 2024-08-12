% -------------------------------------------------------------------------
% function to plot general controller fit, which will then be usable by new
% code that tries to synthesize controller fitting routines for different
% DOFs
% -------------------------------------------------------------------------
function [h, ax] = plotControllerFit(controller_fit_struct)
% -------------------------------
%% params
PI_color = [70,130,180]/255 ; % color of controller fit line
I_color = [0,0,0] ;
P_color = 0.6*[1 1 1] ;
data_color = [0,0,0] ;
pulse_color = [255 238 170]/255 ;

pulse_alpha = 0.4 ;

I_lineStyle = '--' ;
P_lineStyle = '-' ;

lineWidthThick = 2.5 ;
lineWidthThin = 1.5 ;

% -------------------------------
%% read info from control struct
% many the relevant kinematic variables have pertType-specific field names
% NB: this is the NUMBER representation
if isfield(controller_fit_struct,'pertType')
    pertType = controller_fit_struct.pertType ;
elseif isfield(controller_fit_struct,'pitchType')
    pertType = controller_fit_struct.pitchType ;
else
    fprintf('Error: cannot determine pert type -- quitting \n')
    h = [] ;
    ax = [] ;
    return
end

% ----------------------------------------------
% switch how we read struct based on pert type
switch pertType
    case {-1,1}
        % this is a PITCH perturbation
        c_bodyAngle = controller_fit_struct.c_pitch ;
        wingAngleTimes = controller_fit_struct.fwdFlipTimes ;
        wingAngleVals = controller_fit_struct.deltaPhiFront ;
        
        pertTypeStr = 'Pitch' ;
        wingAngleName = '\Delta\phi_{fwd}' ;
        
    case {-2,2}
        % this is a ROLL perturbation
        c_bodyAngle = controller_fit_struct.c_roll ;
        wingAngleTimes = controller_fit_struct.phiAmpTimes ;
        wingAngleVals = controller_fit_struct.phiAmpDiff ;
        
        pertTypeStr = 'Roll' ;
        wingAngleName = '\Delta\Phi_{LR}' ;
        
    case {-3,3}
        % this is a YAW perturbation
        c_bodyAngle = controller_fit_struct.c_yaw ;
        wingAngleTimes = controller_fit_struct.midWingBeatTimes ;
        wingAngleVals = controller_fit_struct.deltaAlpha ;
        
        pertTypeStr = 'Roll' ;
        wingAngleName = '\Delta\Phi_{LR}' ;
        
    otherwise
        fprintf('Error: cannot determine pert type -- quitting \n')
        h = [] ;
        ax = [] ;
        return
end

% also read out common terms
K_i = controller_fit_struct.K_i ;
K_p = controller_fit_struct.K_p ;
deltaT = controller_fit_struct.deltaT ;

% try to get pulse info
if isfield(controller_fit_struct,'pulseTiming')
    % if we know pulse info, read it
    pulseTiming = controller_fit_struct.pulseTiming ;
else
    % otherwise just have to guess...
    pulseTiming = (1e-3).*[0, 7] ;
end

% if there's a constant term in controller fit, load it
if isfield(controller_fit_struct,'K')
    K = controller_fit_struct.K ;
else
    % otherwise set to zero
    K = 0 ;
end
% -------------------------------------------------------------------
%% get derivatives, time range, etc
% get t limits (round to nearest 5 ms)
t_start = 0.005*floor(wingAngleTimes(1)*200) ;
t_end = 0.005*ceil(wingAngleTimes(end)*200) ;
t = linspace(t_start, t_end, 100) ;

% evaluate interpolant for body angle at these points
bodyAngle = c_bodyAngle(t - deltaT) ;
if ismember(pertTypeStr, {'Pitch','Yaw'})
    bodyAngleInit = c_bodyAngle(pulseTiming(1)) ;
    deltaBodyAngle = bodyAngle - bodyAngleInit ;
else
    deltaBodyAngle = bodyAngle ;
end

% also get angular velocity
bodyAngleVel = differentiate(c_bodyAngle, t - deltaT) ;

% combine to get controller fit prediction
controlPred = (K_i * deltaBodyAngle) + (K_p * bodyAngleVel) + K ;

% fit confidence interval
if ismember(pertTypeStr, {'Pitch','Roll'})
    [upper, lower, t_CI] = get_controllerFit_CI(controller_fit_struct) ;
    CIFlag = true ; 
else
    upper = [] ; 
    lower = [] ;
    CIFlag = false ; 
end
% -------------------------------------------------------------------
%% make plot
% figure out some of the limits
ylim = [min(wingAngleVals)-5 , max(wingAngleVals)+5] ;
xlim = 1000.*[t(1), t(end)] ;
tsfvec = 1000.*[pulseTiming(1), pulseTiming(2), pulseTiming(2), ...
    pulseTiming(1), pulseTiming(1)] ; % need to convert to ms
tms = 1000.*t ;  % time in milliseconds

% make figure window
h = figure('PaperPositionMode','auto','Position', [500 500 420 140]) ;
ax = gca ;
hold on

% draw patch for magnetic pulse
avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
hf = fill(ax, tsfvec , avec,'y') ;
set(hf,'facecolor',pulse_color,'facealpha',pulse_alpha,'edgecolor','none') ;
set(hf,'HandleVisibility','off')

% draw controller terms
plot(ax, tms, K_p * bodyAngleVel,...
    'Color',P_color,...
    'LineWidth',lineWidthThin,...
    'LineStyle', P_lineStyle)
plot(ax, tms, K_i * deltaBodyAngle,...
    'Color',I_color,...
    'LineWidth',lineWidthThin,...
    'LineStyle', I_lineStyle)

% draw full controller fit (and CI?)
if CIFlag
    h_CI = fill(ax, 1000.*[t_CI,fliplr(t_CI)],[lower',fliplr(upper')],...
        PI_color,'linestyle','none');
    h_CI.FaceAlpha = pulse_alpha ;
end
plot(ax, tms, controlPred,...
    'Color',PI_color,...
    'LineWidth',lineWidthThick,...
    'LineStyle', '-')

% draw data
plot(ax, 1000*wingAngleTimes, wingAngleVals,'o',...
    'markerfacecolor',data_color,...
    'markeredgecolor',data_color) ;

% ---------------------------------------------
% set axis properties
set(ax,'xlim',xlim,'ylim',ylim)
xlabel('time (ms)')
ylabel(sprintf('%s (deg)', wingAngleName))

end