% -------------------------------------------------------------------------
% function to fit PI controller model to L/R difference in stroke amplitude
%  for roll perturbations NEED TO FILL THIS IN MORE
% -------------------------------------------------------------------------
function controller_fit_struct = fitRollController_dt(rootPath, MovNum, ...
    saveFlag, plotFlag1, plotFlag2, debugFlag1, debugFlag2, smoothFlag,...
    verboseFlag)
% -------------------------------------
%% inputs 
if ~exist('rootPath','var') || isempty(rootPath)
    rootPath = 'D:\Fly Data\VNC Motor Lines\76_15062019\' ; %'D:\Fly Data\VNC Motor Lines\52_19022019\' ;
end
if ~exist('MovNum','var') || isempty(MovNum) 
    MovNum = [] ; %5 ; % 14, 17
end
if ~exist('saveFlag','var') || isempty(saveFlag) 
    saveFlag = false ;
end
if ~exist('plotFlag1','var') || isempty(plotFlag1) 
    plotFlag1 = false ;
end
if ~exist('plotFlag2','var') || isempty(plotFlag2) 
    plotFlag2 = true ;
end
if ~exist('debugFlag1','var') || isempty(debugFlag1) 
    debugFlag1 = true ;
end
if ~exist('debugFlag2','var') || isempty(debugFlag2) 
    debugFlag2 = true ;
end
if ~exist('smoothFlag','var') || isempty(smoothFlag) 
    smoothFlag = true ;
end
if ~exist('verboseFlag','var') || isempty(verboseFlag) 
    verboseFlag = true ;
end

% -------------------------------------------
%% params
fit_params = setControllerFitParams() ; 
% initial guess
K_i_guess = fit_params.K_i_guessRoll ;%0.3 ; %0.3 ; %unitless
K_p_guess = fit_params.K_p_guessRoll ; %seconds
gainGuess = [K_i_guess, K_p_guess] ;  %K_i, K_p, K (from paper)

% space of time delays to search 
deltaT_min = fit_params.deltaT_min ; %0.003
deltaT_max = fit_params.deltaT_max ; %0.013
deltaT_N = fit_params.deltaT_N ; %60
deltaTGuess = fit_params.deltaTGuessRoll ;
deltaT_space = linspace(deltaT_min, deltaT_max, deltaT_N) ;

% range over which to fit
fitRangeMS = fit_params.rollFitRangeMS ; 

% initialize controller fit structure
controller_fit_struct = struct() ; 

% -----------------------------------------------------------------
%% get data
%full data structure
[data, dataPath, pertType, ExprNum] = loadPertDataStruct(rootPath, MovNum) ;
data.manualCorrRangeMS = fitRangeMS ;

%data to fit controller model to
[phiAmpDiff, c_roll, phiAmpTimes] = getRollControllerData(data, debugFlag1) ;

if smoothFlag
    phiAmpDiff = smooth(phiAmpDiff,3)' ; %smooth(phiAmpDiff,3)' ; 
end
%phiAmpDiff =phiAmpDiff + 10 ;  
N = length(phiAmpTimes) ;

% -----------------------------------------------------------------
%% loop through deltaT and fit gains
resNormArray = nan(deltaT_N,1) ;
residualArray = nan(deltaT_N,N) ;
K_i_Array = nan(deltaT_N,1) ;
K_p_Array = nan(deltaT_N,1) ;
J_cell = cell(deltaT_N,1) ;
Rsq_array = nan(deltaT_N,1) ;

%cc = 1 ;
for k = 1:deltaT_N
    deltaT = deltaT_space(k) ;
    [K_i, K_p,resnorm, residuals, J, Rsq_adj] = ...
        fitRollControllerGains(phiAmpTimes, phiAmpDiff, deltaT, c_roll,...
        gainGuess, plotFlag1)  ;
    if plotFlag1
        title(['\DeltaT = ' num2str(deltaT)])
    end
    
    resNormArray(k) = resnorm ;
    residualArray(k,:) = residuals ;
    K_i_Array(k) = K_i ;
    K_p_Array(k) = K_p ;
    J_cell{k} = J ;
    Rsq_array(k) = Rsq_adj ; 
    %cc = cc + 1 ;
end

% -----------------------------------------------------------------
%% estimate jacobian for deltaT
J_deltaT_Array = zeros(N, deltaT_N) ;

for i = 1:N
    residuals = residualArray(:,i) ;
    c_residual = fit(deltaT_space', residuals, 'cubicinterp') ;
    
    residual_derivative = differentiate(c_residual, deltaT_space) ;
    J_deltaT_Array(i,:) = residual_derivative' ;
end

% ------------------------------------------------------
%% find best fit
% previously was just: [~, bestFitInd] = min(resNormArray) ;
[rms_minima, rms_minima_ind] = findpeaks(-1*resNormArray) ; 
if isempty(rms_minima)
    [~, bestFitInd] = min(resNormArray) ;
elseif length(rms_minima) == 1 
    bestFitInd = rms_minima_ind ; 
elseif length(rms_minima) > 1
    rms_minima_t = deltaT_space(rms_minima_ind) ; 
    [~, temp_idx] = min(abs(rms_minima_t - deltaTGuess)) ;
    bestFitInd = rms_minima_ind(temp_idx) ;
end

K_i = K_i_Array(bestFitInd) ;
K_p = K_p_Array(bestFitInd) ;
deltaT = deltaT_space(bestFitInd) ;
resnorm = resNormArray(bestFitInd) ;
rmse = sqrt(resnorm/N) ;
Rsq_adj = Rsq_array(bestFitInd) ; 

% --------------------------------------------------------
%% use jacobian to estimate confidence intervals
Jacobian_gains = J_cell{bestFitInd} ;
Jacobian_gains = full(Jacobian_gains) ;

J_deltaT = J_deltaT_Array(:, bestFitInd) ;
Jacobian = [Jacobian_gains, J_deltaT] ;

varp = resnorm*inv(Jacobian'*Jacobian)/N ;
stdp = sqrt(diag(varp));
%stdp = 100*stdp'./p;

% ------------------------------------------------------
%% diagnostic plots?
if debugFlag2
    figure ;
    hold on
    stem(deltaT_space, resNormArray)
    plot(deltaT_space(bestFitInd), resNormArray(bestFitInd), 'rx')
end

if plotFlag2
    [~,~, ~, ~, ~] = ...
        fitRollControllerGains(phiAmpTimes, phiAmpDiff, ...
        deltaT_space(bestFitInd), c_roll, gainGuess, true)  ;
end

% -------------------------------------------------------
%% print results?
if verboseFlag
    fprintf('K_i = %f +/- %f \n', K_i, stdp(1))
    fprintf('K_p = %f +/- %f \n', K_p, stdp(2))
    fprintf('DT = %f +/- %f \n', deltaT, stdp(3))
   
    fprintf('RMSE: \n %f \n', rmse)
    fprintf('R squared (adjusted): \n %f \n', Rsq_adj)
end

% ----------------------------------------------------------
%% assign confidence intervals
K_i_CI = stdp(1) ; 
K_p_CI = stdp(2) ; 
deltaT_CI = stdp(3) ;

% ------------------------------------------------------------
%% save to struct
controller_fit_struct.ExprNum = ExprNum ;
controller_fit_struct.MovNum = MovNum ;
controller_fit_struct.K_i = K_i ;
controller_fit_struct.K_p = K_p ;
controller_fit_struct.deltaT = deltaT ;
controller_fit_struct.K_i_CI = K_i_CI ;
controller_fit_struct.K_p_CI = K_p_CI ;
controller_fit_struct.deltaT_CI = deltaT_CI ;
controller_fit_struct.rms = rmse ;
controller_fit_struct.Rsq_adj = Rsq_adj ;
controller_fit_struct.c_roll = c_roll ;
%controller_fit_struct.sp_phiR = sp_phiR ;
%controller_fit_struct.sp_phiL = sp_phiL ;
%controller_fit_struct.pitchEstErr = pitchEstErr ;
%controller_fit_struct.phiEstErr = phiEstErr ;
controller_fit_struct.phiAmpTimes = phiAmpTimes ;
%controller_fit_struct.backFlipTimes = backFlipTimes ;
controller_fit_struct.phiAmpDiff = phiAmpDiff ;
controller_fit_struct.pertType = pertType ;
controller_fit_struct.rootPath = rootPath ;

if saveFlag
    save(fullfile(dataPath, 'controller_fit_struct_LM.mat'),...
        'controller_fit_struct')
end

end
