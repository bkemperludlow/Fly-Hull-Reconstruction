% -------------------------------------------------------------------------
% function to perform linear least squares fitting for controller gains at
% a given vaule of deltaT for arbitrary controller (i.e. angular DOF) fit
% -------------------------------------------------------------------------
function [x, resnorm, residual, J, Rsq_adjusted] = ...
    fitControllerGains(wingAngleTimes, wingAngleVals, deltaT, ...
    c_bodyAngle, pertType, paramGuess, plotFlag, pulseTiming, KFlag) 
% ---------------------------------
%% inputs
if ~exist('pulseTiming','var') || isempty(pulseTiming)
    % two element vector giving pulse start and stop in seconds
    pulseTiming = [0.0, 0.007] ; 
end
if ~exist('KFlag','var') || isempty(KFlag)
    % fit for a constant "K" in the controller equation?
    KFlag = false; 
end

% ------------------------------------------------
%% prepare kinematic data
bodyAngle = c_bodyAngle(wingAngleTimes - deltaT) ; 

% if we're doing a pitch or yaw controller fit, should subtract off initial
% value of body angle
if ismember(pertType, {'Pitch', 'Yaw'})
    initial_t =  pulseTiming(1) ; % 0 ; % (-0.0025 : 1.25e-4 : 0) ; %
    bodyAngleInit = mean(c_bodyAngle(initial_t)); % c_pitch(0) ; %
    deltaBodyAngle = bodyAngle - bodyAngleInit ;
else
    deltaBodyAngle = bodyAngle ; 
end

% also get body angle velocity
bodyAngleVelocity = differentiate(c_bodyAngle, wingAngleTimes - deltaT) ;

% combine angular disp and velocity into one array
bodyData = [deltaBodyAngle' ; bodyAngleVelocity'] ; 

% -------------------------------------------------------------------
%% define controller equations and set solver options
% define controller equation (include constant or not -- probably NOT)
if KFlag
    % this controller includes a constant term (x(3))
    controllerEq = @(x, xdata) x(1)*xdata(1,:) + x(2)*xdata(2,:) + x(3);
    
    % if we are fitting for a constant, make sure param guess includes an
    % initial estimate for K
    if length(paramGuess) < 3
       paramGuess(end+1) = 0 ;  
    end
else
    % this controller equation just contains 2 gain coefficients (K_i, K_p)
    controllerEq = @(x, xdata) x(1)*xdata(1,:) + x(2)*xdata(2,:); 
end

% set solver options
options = optimoptions('lsqcurvefit',...
    'Algorithm','levenberg-marquardt',...
    'display','off',...
    'FiniteDifferenceType','central',...
    'FunctionTolerance',1e-10) ; %, ...
%     'FiniteDifferenceStepSize', 1e-5); 
%'trust-region-reflective' or 'levenberg-marquardt'

% define bounds (can't use bounds with 'levenberg-marquardt')
lb = [] ;%[-2 -0.1 -4]; %
ub = [] ; %[ 2  0.1  4];

% ---------------------------------------------------------------------
%% perform fit and assign outputs 
[x, resnorm, residual, ~, ~, ~, J] = ...
    lsqcurvefit(controllerEq, paramGuess, bodyData,wingAngleVals',...
    lb,ub,options) ;

% ----------------------------------------------------------------
%% calculate adjusted coefficient of determination (adjusted R^2)
SS_tot = sum((wingAngleVals - mean(wingAngleVals)).^2) ; 
SS_reg = resnorm ; 
n_obs = length(wingAngleVals) ; 
n_params = length(x) ; 
Rsq_adjusted = 1 - ((n_obs - 1)/(n_obs - n_params))*(SS_reg/SS_tot) ; 

% ----------------------------------------------------------------------
%% plot results?
if plotFlag
    % make temporary controller fit struct and use that with extant
    % function
    fprintf('Plotting in gain fit function is under construction! \n')
    %{
    t = linspace(wingAngleTimes(1), wingAngleTimes(end), 100) ; 
    deltaBodyAngle_cont = c_bodyAngle(t - deltaT) - bodyAngleInit ; %continuous
    bodyAngleVelocity_cont = differentiate(c_bodyAngle, t - deltaT) ;
    controlPred = K_i * deltaBodyAngle_cont + K_p * bodyAngleVelocity_cont ; %+ K ; 
    ylim = [min(wingAngleVals)-5 , max(wingAngleVals)+5] ;
    tsfvec = 1000.*[pulseTiming(1), pulseTiming(2), pulseTiming(2), ...
        pulseTiming(1), pulseTiming(1)] ; % need to convert to ms
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
    plot(1000*t,K_p * bodyAngleVelocity_cont,'Color',0.6*[1 1 1],'LineWidth',1.5)
    plot(1000*t,K_i * deltaBodyAngle_cont,'k--','LineWidth',1.5)
    errorbar(1000*wingAngleTimes,wingAngleVals, 2*ones(size(wingAngleVals)),'ko','markerfacecolor','k') ;
    %plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    axis tight ;
    set(gca, 'xlim', 1000*[t(1), t(end)])
    set(gca,'ylim',ylim)
    xlabel('Time [ms]')
    ylabel('\Delta wing angle [deg]')
    %legend({'PI','P','I','Data'})
   %}

end