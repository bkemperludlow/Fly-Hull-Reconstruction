% -------------------------------------------------------------------------
% function to perform linear least squares fitting for controller gains at
% a given vaule of deltaT for pitch controller fits (NO K!)
% -------------------------------------------------------------------------
function [K_i, K_p, resnorm, residual, J, Rsq_adjusted] = ...
    fitPitchControllerGains_noK(fwdFlipTimes, deltaPhiFront, deltaT, ...
    c_pitch, paramGuess, plotFlag) 
% ------------------------------------------------
%% prepare kinematic data

bodyPitch = c_pitch(fwdFlipTimes - deltaT) ; 
initial_t =  0 ; % (-0.0025 : 1.25e-4 : 0) ; %
theta_0 = mean(c_pitch(initial_t)); % c_pitch(0) ; %  
pitchVelocity = differentiate(c_pitch, fwdFlipTimes - deltaT) ;
deltaBodyPitch = bodyPitch - theta_0 ; 
bodyData = [deltaBodyPitch' ; pitchVelocity'] ; 

% -------------------------------------------------------------------
%% define controller equations and set solver options
controllerEq = @(x, xdata) x(1)*xdata(1,:) + x(2)*xdata(2,:); 

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',...
    'display','off','FiniteDifferenceType','central',...
    'FunctionTolerance',1e-10); 
%'trust-region-reflective' or 'levenberg-marquardt'
lb = [] ;%[-2 -0.1 -4]; %
ub = [] ; %[ 2  0.1  4];

% ---------------------------------------------------------------------
%% perform fit and assign outputs 
[x, resnorm, residual, ~, ~, ~, J] = ...
    lsqcurvefit(controllerEq, paramGuess, bodyData,deltaPhiFront,...
    lb,ub,options) ;

K_i = x(1) ;
K_p = x(2) ;

% ----------------------------------------------------------------
%% calculate adjusted coefficient of determination (adjusted R^2)
SS_tot = sum((deltaPhiFront - mean(deltaPhiFront)).^2) ; 
SS_reg = resnorm ; 
n_obs = length(deltaPhiFront) ; 
n_params = length(x) ; 
Rsq_adjusted = 1 - ((n_obs - 1)/(n_obs - n_params))*(SS_reg/SS_tot) ; 

% ----------------------------------------------------------------------
%% plot results?
if plotFlag
    t = linspace(fwdFlipTimes(1), fwdFlipTimes(end), 100) ; 
    deltaBodyPitch_cont = c_pitch(t - deltaT) - theta_0 ; %continuous
    pitchVelocity_cont = differentiate(c_pitch, t - deltaT) ;
    controlPred = K_i * deltaBodyPitch_cont + K_p * pitchVelocity_cont ; %+ K ; 
    ylim = [min(deltaPhiFront)-5 , max(deltaPhiFront)+5] ;
    tsfvec = [0 7 7 0 0] ;
    %patchColor = [1 1 1 ] * 0.8;
    plotColor = [70,130,180]/255 ; 
        
    figure ; 
    hold on
       
    set(gcf, 'Position', [500 500 420 140]);
    set(gcf,'PaperPositionMode','auto')
    
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',1,'edgecolor','none') ;
    set(hf,'HandleVisibility','off')
    
    plot(1000*t,controlPred,'Color',plotColor,'LineWidth',2.5)
    plot(1000*t,K_p * pitchVelocity_cont,'Color',0.6*[1 1 1],'LineWidth',1.5)
    plot(1000*t,K_i * deltaBodyPitch_cont,'k--','LineWidth',1.5)
    errorbar(1000*fwdFlipTimes,deltaPhiFront, 2*ones(size(deltaPhiFront)),'ko','markerfacecolor','k') ;
    %plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    axis tight ;
    set(gca, 'xlim', 1000*[t(1), t(end)])
    set(gca,'ylim',ylim)
    xlabel('Time [ms]')
    ylabel('\Delta \Phi_{front} [deg]')
    %legend({'PI','P','I','Data'})
   

end