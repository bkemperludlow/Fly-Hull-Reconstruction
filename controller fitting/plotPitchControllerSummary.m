% -------------------------------------------------------------------------
% function to plot summary of controller fit, including body pitch, pitch
% velocity, fwd stroke angle, fit and parameters
% -------------------------------------------------------------------------
function h_main = plotPitchControllerSummary(controller_fit_struct)
% -----------------------
%% params
plotColor = [70,130,180]/255 ; 
patchColor = [255 238 170]/255;
figPosition = [2341, 86, 453, 780] ; 
% duration of magnetic pulse
pulseDurationMS = 7 ; % ms

% -------------------------
%% read in relevant data
if isfield(controller_fit_struct,'K')
    KFlag = true ; 
    K = controller_fit_struct.K ; 
    K_CI = controller_fit_struct.K_CI ; 
else
    KFlag = false ; 
    K = 0 ; 
    K_CI = 0 ; 
end

% kinematics
fwdFlipTimes = controller_fit_struct.fwdFlipTimes ; 
deltaPhiFront = controller_fit_struct.deltaPhiFront ;
c_pitch = controller_fit_struct.c_pitch ;

% movie info
if isfield(controller_fit_struct,'driver') && ...
        isfield(controller_fit_struct,'effector')
    driver = controller_fit_struct.driver ;
    effector = controller_fit_struct.effector ;
else
    driver = '?' ;
    effector = '?' ;
end
ExprNum = controller_fit_struct.ExprNum ;
MovNum = controller_fit_struct.MovNum ;
movieStr = ['Expr_' num2str(ExprNum) '_mov_' num2str(MovNum,'%03d')] ;

% fit params
K_i = controller_fit_struct.K_i ; 
K_p = controller_fit_struct.K_p ; 
deltaT = controller_fit_struct.deltaT ; 

K_i_CI = controller_fit_struct.K_i_CI ; 
K_p_CI = controller_fit_struct.K_p_CI ; 
deltaT_CI = controller_fit_struct.deltaT_CI ; 
rmse = controller_fit_struct.rms ; 

% ---------------------------------------------------------------------
%% curves to plot 
% time
t = linspace(fwdFlipTimes(1), fwdFlipTimes(end), 100) ; 
tms = 1000*t ;  
xlim = [tms(1), tms(end)] ; 

% body pitch
theta_0 = c_pitch(0) ;
bodyPitch = c_pitch(t)  ;
deltaBodyPitch = bodyPitch - theta_0 ; 

% pitch velocity
pitchVelocity = differentiate(c_pitch, t) ;
pitchVelocity_dt = differentiate(c_pitch, t - deltaT) ; 

% controller terms
I_term = K_i.*(c_pitch(t - deltaT) - theta_0) ; 
P_term = K_p.*(pitchVelocity_dt) ; 

PI_model = I_term + P_term + K ; 

% ---------------------------------------------------------------------
%% create plot
h_main = figure('PaperPositionMode','auto','Position',figPosition) ; 
% yyaxis color order
left_color = [0, 0, 1] ; 
right_color = [0, 0, 0] ; 
set(h_main,'defaultAxesColorOrder',[left_color; right_color]);

% ----------------------
% BODY PITCH
% ----------------------
ax_pitch = subplot(4,1,1) ; 
hold on
yyaxis left

% pert patch
tsfvec = [0 pulseDurationMS pulseDurationMS 0 0] ;
ylim_curr = [min(deltaBodyPitch) - 5, max(deltaBodyPitch) + 5] ; 
avec = [ylim_curr(1) ylim_curr(1) ylim_curr(2) ylim_curr(2) ylim_curr(1)] ;
hf = fill(tsfvec , avec,'y') ;
set(hf,'facecolor',patchColor,'facealpha',0.5,'edgecolor','none') ;
set(hf,'HandleVisibility','off')

% left axis
plot(tms, deltaBodyPitch, 'b-','LineWidth',1.5)
ylabel('\Delta \theta')

yyaxis right
plot(1000*(fwdFlipTimes-deltaT), deltaPhiFront, 'o--', ...
    'Color', 0.6*[1,1,1], 'MarkerFaceColor',0.6*[1,1,1])
plot(1000*fwdFlipTimes, deltaPhiFront, 'ko-','MarkerFaceColor','k')
ylabel('\Delta \Phi_{front}')

set(ax_pitch, 'xlim', xlim)
title(movieStr,'Interpreter', 'none')

% ----------------------
% PITCH VELOCITY
% ----------------------
ax_vel = subplot(4,1,2) ; 
hold on

% pert patch
tsfvec = [0 pulseDurationMS pulseDurationMS 0 0] ;
ylim_curr = [min(pitchVelocity) - 500, max(pitchVelocity) + 500] ; 
avec = [ylim_curr(1) ylim_curr(1) ylim_curr(2) ylim_curr(2) ylim_curr(1)] ;
hf = fill(tsfvec , avec,'y') ;
set(hf,'facecolor',patchColor,'facealpha',0.5,'edgecolor','none') ;
set(hf,'HandleVisibility','off')

% left axis
yyaxis left
plot(tms, pitchVelocity, 'b-','LineWidth',1.5)
ylabel('Pitch Velocity')

yyaxis right
hold on
plot(1000*(fwdFlipTimes-deltaT), deltaPhiFront, 'o--', ...
    'Color', 0.6*[1,1,1], 'MarkerFaceColor',0.6*[1,1,1])
plot(1000*fwdFlipTimes, deltaPhiFront, 'ko-','MarkerFaceColor','k')
ylabel('\Delta \Phi_{front}')

set(ax_vel, 'xlim', xlim)

% ----------------------
% CONTROLLER FIT
% ----------------------
ax_fit = subplot(4,1,3) ; 
hold on 

% pert patch
tsfvec = [0 pulseDurationMS pulseDurationMS 0 0] ;
ylim_curr = [min(deltaPhiFront) - 5, max(deltaPhiFront) + 5] ; 
avec = [ylim_curr(1) ylim_curr(1) ylim_curr(2) ylim_curr(2) ylim_curr(1)] ;
hf = fill(tsfvec , avec,'y') ;
set(hf,'facecolor',patchColor,'facealpha',0.5,'edgecolor','none') ;
set(hf,'HandleVisibility','off')

% controller terms
plot(tms,P_term,'Color',0.6*[1 1 1],'LineWidth',1.0)
plot(tms,I_term,'k--','LineWidth',1.0)
plot(tms,PI_model,'Color',plotColor,'LineWidth',1.5)
errorbar(1000*fwdFlipTimes,deltaPhiFront, 2*ones(size(deltaPhiFront)),'ko','markerfacecolor','k') ;
%plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
axis tight ;
set(ax_fit, 'xlim', xlim)
xlabel('Time (ms)')
ylabel('\Delta \Phi_{front} (deg)')

% ----------------------
% TEXT SUMMARY OF FIT
% ----------------------
ax_text = subplot(4,1,4) ;
txt_gtype = text(ax_text, 0.1 , 0.85, [effector ' > ' driver]) ;
txt_ki = text(ax_text, 0.1, 0.6, ...
    ['K_i = ' num2str(K_i) ' +/- ' num2str(K_i_CI)]) ;
txt_kp = text(ax_text, 0.1, 0.40, ...
    ['K_p = ' num2str(K_p) ' +/- ' num2str(K_p_CI)]) ;
txt_dt = text(ax_text, 0.1, 0.20, ...
    ['\Delta T = ' num2str(deltaT) ' +/- ' num2str(deltaT_CI)]) ;
if KFlag
    txt_k = text(ax_text, 0.1, 0.0, ...
        ['K = ' num2str(K) ' +/- ' num2str(K_CI)]) ;
end
txt_rmse = text(ax_text, 0.1, -0.10, ['RMSE = ' num2str(rmse)]) ;
if isfield(controller_fit_struct,'Rsq_adj')
    Rsq_adj = controller_fit_struct.Rsq_adj ;
    txt_rsq = text(ax_text, 0.1, -0.25, ['R^2 (adj) = ' num2str(Rsq_adj)]) ;
end
axis off



end