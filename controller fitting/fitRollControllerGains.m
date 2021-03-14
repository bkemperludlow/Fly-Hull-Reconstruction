% -------------------------------------------------------------------------
% function to perform linear least squares fitting for controller gains at
% a given vaule of deltaT for ROLL controller fits 
% -------------------------------------------------------------------------
function [K_i, K_p, resnorm, residual, J, Rsq_adjusted] = ...
    fitRollControllerGains(phiAmpTimes, phiAmpDiff, deltaT, c_roll, ...
    paramGuess, plotFlag) 
% --------------------------------------------------------
%% prepare kinematic data
bodyRoll = c_roll(phiAmpTimes - deltaT) ; 
rollVelocity = differentiate(c_roll, phiAmpTimes - deltaT) ;
bodyData = [bodyRoll' ; rollVelocity'] ; 

% --------------------------------------------------------------
%% define controller equation and solver options
controllerEq = @(x, xdata) x(1)*xdata(1,:) + x(2)*xdata(2,:) ; 

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',...
    'display','off','FiniteDifferenceType','central','FunctionTolerance',1e-10); 
%'trust-region-reflective' or 'levenberg-marquardt'
lb = [] ;%[-2 -0.1 -4]; %
ub = [] ; %[ 2  0.1  4];

% --------------------------------------------------------------
%% run solver and assign outputs
[x, resnorm, residual, ~, ~, ~, J] = ...
    lsqcurvefit(controllerEq, paramGuess, bodyData,phiAmpDiff,lb,ub,options) ;

K_i = x(1) ;
K_p = x(2) ;

% ----------------------------------------------------------------
%% calculate adjusted coefficient of determination (adjusted R^2)
SS_tot = sum((phiAmpDiff - mean(phiAmpDiff)).^2) ; 
SS_reg = resnorm ; 
n_obs = length(phiAmpDiff) ; 
n_params = length(x) ; 
Rsq_adjusted = 1 - ((n_obs - 1)/(n_obs - n_params))*(SS_reg/SS_tot) ; 

% -----------------------------------------------------------------
%% plot results?
if plotFlag
    t = linspace(phiAmpTimes(1), phiAmpTimes(end), 100) ; 
    bodyRoll_cont = c_roll(t - deltaT) ; %continuous
    rollVelocity_cont = differentiate(c_roll, t - deltaT) ;
    controlPred = K_i * bodyRoll_cont + K_p * rollVelocity_cont ; 
    ylim = [min(phiAmpDiff)-5 , max(phiAmpDiff)+5] ;
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
    plot(1000*t,K_p * rollVelocity_cont,'Color',0.6*[1 1 1],'LineWidth',1.5)
    plot(1000*t,K_i * bodyRoll_cont,'k--','LineWidth',1.5)
    errorbar(1000*phiAmpTimes,phiAmpDiff, 2*ones(size(phiAmpDiff)),'ko','markerfacecolor','k') ;
    %plotWingstrokeBackground(gca, backFlipTimes*1000, phiAmpTimes*1000, patchColor, true);
    axis tight ;
    set(gca, 'xlim', 1000*[t(1), t(end)])
    set(gca,'ylim',ylim)
    xlabel('Time [ms]')
    ylabel('\Delta \Phi_{LR} [deg]')
    %legend({'PI','P','I','Data'})
   

end