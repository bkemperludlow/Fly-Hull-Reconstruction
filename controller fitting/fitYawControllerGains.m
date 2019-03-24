function [K_i, K_p, resnorm, residual, J] = ...
    fitYawControllerGains(midWingBeatTimes, deltaAlphaMean, deltaT, c_yaw, paramGuess, plotFlag) 

%x = [ K_i, K_p, K] 

bodyYaw = c_yaw(midWingBeatTimes - deltaT) ; 
phib_0 = c_yaw(0) ;
yawVelocity = differentiate(c_yaw, midWingBeatTimes - deltaT) ;
deltaBodyYaw = bodyYaw - phib_0 ; 
bodyData = [deltaBodyYaw' ; yawVelocity'] ; 

controllerEq = @(x, xdata) x(1)*xdata(1,:) + x(2)*xdata(2,:) ; 

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt',...
    'display','off','FiniteDifferenceType','central','FunctionTolerance',1e-10); 
%'trust-region-reflective' or 'levenberg-marquardt'
lb = [] ;%[-2 -0.1 -4]; %
ub = [] ; %[ 2  0.1  4];

% options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
%     'display','off','FiniteDifferenceType','central','FunctionTolerance',1e-10); 
% %'trust-region-reflective' or 'levenberg-marquardt'
% lb = [-2 0 -4]; %
% ub =[ 2 0.1  4];

[x, resnorm, residual, ~, ~, ~, J] = ...
    lsqcurvefit(controllerEq, paramGuess, bodyData,deltaAlphaMean',lb,ub,options) ;

K_i = x(1) ;
K_p = x(2) ;


if plotFlag
    t = linspace(midWingBeatTimes(1), midWingBeatTimes(end), 100) ; 
    deltaBodyYaw_cont = c_yaw(t - deltaT) - phib_0 ; %continuous
    yawVelocity_cont = differentiate(c_yaw, t - deltaT) ;
    controlPred = K_i * deltaBodyYaw_cont + K_p * yawVelocity_cont  ; 
    ylim = [min(deltaAlphaMean)-5 , max(deltaAlphaMean)+5] ;
    tsfvec = [0 7 7 0 0] ;
    %patchColor = [1 1 1 ] * 0.8;
    plotColor = [0,109,44]/255 ; 
        
    figure ; 
    hold on
       
    set(gcf, 'Position', [500 500 420 140]);
    set(gcf,'PaperPositionMode','auto')
    
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',1,'edgecolor','none') ;
    set(hf,'HandleVisibility','off')
    
    plot(1000*t,controlPred,'Color',plotColor,'LineWidth',2.5)
    plot(1000*t,K_p * yawVelocity_cont,'Color',0.6*[1 1 1],'LineWidth',1.5)
    plot(1000*t,K_i * deltaBodyYaw_cont,'k--','LineWidth',1.5)
    errorbar(1000*midWingBeatTimes,deltaAlphaMean, 2*ones(size(deltaAlphaMean)),'ko','markerfacecolor','k') ;
    %plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    axis tight ;
    set(gca, 'xlim', 1000*[t(1), t(end)])
    set(gca,'ylim',ylim)
    xlabel('Time [ms]')
    ylabel('\Delta \alpha_{LR} [deg]')
    %legend({'PI','P','I','Data'})
   

end