% -------------------------------------------------------------------------
% function to fit PI controller model to front stroke angle for pitch
% perturbations. NEED TO FILL THIS IN MORE
% -------------------------------------------------------------------------
function controller_fit_struct = fitPitchController_dt(rootPath, MovNum, ...
    saveFlag, KFlag, plotFlag1, plotFlag2, debugFlag1, debugFlag2, ...
    smoothFlag, verboseFlag, pulseDurationIn)
% -------------------------------------
%% inputs
if ~exist('rootPath','var') || isempty(rootPath)
    rootPath = 'D:\Fly Data\VNC Motor Lines\96_16102019\'   ; %'D:\Fly Data\VNC Motor Lines\52_19022019\' ;
end
if ~exist('MovNum','var') || isempty(MovNum)
    MovNum = 78 ; %5 ; % 14, 17
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = false  ;
end
if ~exist('KFlag','var') || isempty(KFlag)
    KFlag = false ;
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
    smoothFlag = false ;
end
if ~exist('verboseFlag','var') || isempty(verboseFlag)
    verboseFlag = true ;
end
if ~exist('pulseDurationIn','var') || isempty(pulseDurationIn)
    pulseDurationIn = 0.007 ;
end
% ------------------------------------------------
%% params
fit_params = setControllerFitParams(pulseDurationIn) ;
% initial guesses for controller gains
K_i_guess = fit_params.K_i_guessPitch ;
K_p_guess = fit_params.K_p_guessPitch ;
gainGuess = [K_i_guess, K_p_guess] ; %K_i, K_p, K (from paper)
if KFlag
    K_guess = fit_params.K_guessPitch ;
    gainGuess = [gainGuess, K_guess] ;
end
% space of time delays to search
deltaT_min = fit_params.deltaT_min ; %0.003
deltaT_max =  fit_params.deltaT_max ; %0.013
deltaT_N = fit_params.deltaT_N ; %60
deltaTGuess = fit_params.deltaTGuessPitch ;

% time window in which to fit controller
fitRangeMS = fit_params.pitchFitRangeMS ; %  [-10, 45] ; %

% duration of magnetic perturbation
pulseDuration = fit_params.pulseDuration ; % in seconds

% initialize data structure
controller_fit_struct = struct() ;

% ------------------------------------------------------------------------
%% get data
%full data structure
[data, dataPath, pertType, ExprNum] = loadPertDataStruct(rootPath, MovNum) ;
%data.oneWing = 'R' ;
data.manualCorrRangeMS = fitRangeMS ; % need to pass this to 'getPitchControllerData.m' but it doesn't refer to manual correction any longer
data.pulseDuration = pulseDuration ;

%data to fit controller model to
[deltaPhiFront, c_pitch,fwdFlipTimes, backFlipTimes] = ...
    getPitchControllerData(data, debugFlag1) ;
% [deltaPhiFront, c_pitch, c_vel, fwdFlipTimes, backFlipTimes] = ...
%     getPitchControllerData_v3(data, debugFlag1) ;

if smoothFlag
    deltaPhiFront = smooth(deltaPhiFront,3)' ;
end

N = length(fwdFlipTimes) ;

% -------------------------------------------------------------------------
%% loop through deltaT and fit gains

% make discretized t space
deltaT_space = linspace(deltaT_min, deltaT_max, deltaT_N) ;

% initialize arrays for output
resNormArray = nan(deltaT_N,1) ;
residualArray = nan(deltaT_N,N) ;
K_i_Array = nan(deltaT_N,1) ;
K_p_Array = nan(deltaT_N,1) ;
K_Array = nan(deltaT_N,1) ;
J_cell = cell(deltaT_N,1) ;
Rsq_array = nan(deltaT_N, 1) ;

%cc = 1 ;
for k = 1:deltaT_N
    deltaT = deltaT_space(k) ;
    
    if KFlag
        [K_i, K_p, K, resnorm, residuals, J, Rsq_adj] = ...
            fitPitchControllerGains(fwdFlipTimes, deltaPhiFront, deltaT,...
            c_pitch, gainGuess, plotFlag1)  ;
        K_Array(k) = K ;
    else
        [K_i, K_p, resnorm, residuals, J, Rsq_adj] = ...
            fitPitchControllerGains_noK(fwdFlipTimes, deltaPhiFront, deltaT, ...
            c_pitch, gainGuess, plotFlag1) ;
%         [K_i, K_p, resnorm, residuals, J, Rsq_adj] = ...
%             fitPitchControllerGains_noK_v2(fwdFlipTimes, deltaPhiFront, deltaT, ...
%             c_pitch, c_vel, gainGuess, plotFlag1) ;
    end
    
    
    resNormArray(k) = resnorm ;
    residualArray(k,:) = residuals ;
    K_i_Array(k) = K_i ;
    K_p_Array(k) = K_p ;
    Rsq_array(k) = Rsq_adj ;
    J_cell{k} = J ;
    
    if plotFlag1
        title(['\DeltaT = ' num2str(deltaT)])
    end
    %cc = cc + 1 ;
end

% -------------------------------------------------------------------------
%% estimate jacobian for deltaT
J_deltaT_Array = zeros(N, deltaT_N) ;

for i = 1:N
    residuals = residualArray(:,i) ;
    c_residual = fit(deltaT_space', residuals, 'cubicinterp') ;
    
    residual_derivative = differentiate(c_residual, deltaT_space) ;
    J_deltaT_Array(i,:) = residual_derivative' ;
end

% -------------------------------------------------------------------------
%% find best fit
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
if KFlag
    K = K_Array(bestFitInd) ;
end
deltaT = deltaT_space(bestFitInd) ;
resnorm = resNormArray(bestFitInd) ;
rmse = sqrt(resnorm/N) ;

Rsq_adj = Rsq_array(bestFitInd) ; % adjusted coefficient of determination
% --------------------------------------------------
%% use jacobian to estimate confidence intervals
Jacobian_gains = J_cell{bestFitInd} ;
Jacobian_gains = full(Jacobian_gains) ;

J_deltaT = J_deltaT_Array(:, bestFitInd) ;
Jacobian = [Jacobian_gains, J_deltaT] ;

varp = resnorm*inv(Jacobian'*Jacobian)/N ;
stdp = sqrt(diag(varp));
%stdp = 100*stdp'./p;

% -------------------------------------------------------------------------
%% display results
if debugFlag2
    figure ;
    hold on
    stem(deltaT_space, resNormArray)
    plot(deltaT_space(bestFitInd), resNormArray(bestFitInd), 'rx')
end

if plotFlag2
    if KFlag
        [~,~, ~, ~, ~, ~] = ...
            fitPitchControllerGains(fwdFlipTimes, deltaPhiFront,...
            deltaT_space(bestFitInd), c_pitch, gainGuess, true)  ;
    else
        [~, ~, ~, ~, ~] = ...
            fitPitchControllerGains_noK(fwdFlipTimes, deltaPhiFront, deltaT, ...
            c_pitch, gainGuess, true) ;
%         [~, ~, ~, ~, ~] = ...
%             fitPitchControllerGains_noK_v2(fwdFlipTimes, deltaPhiFront, deltaT, ...
%             c_pitch, c_vel, gainGuess, true) ;
    end
    
    
end

% ----------------------------------------------------
%% print results in command window?
if verboseFlag
    fprintf('K_i = %f +/- %f \n', K_i, stdp(1))
    fprintf('K_p = %f +/- %f \n', K_p, stdp(2))
    if KFlag
        fprintf('K = %f +/- %f \n', K, stdp(3))
        fprintf('DT = %f +/- %f \n', deltaT, stdp(4))
    else
        fprintf('DT = %f +/- %f \n', deltaT, stdp(3))
    end
    
    fprintf('RMSE: \n %f \n', rmse)
    fprintf('R squared (adjusted): \n %f \n', Rsq_adj)
end

% -----------------------------------------------------
%% store confidence intervals
if KFlag
    K_i_CI = stdp(1) ;
    K_p_CI = stdp(2) ;
    K_CI = stdp(3) ;
    deltaT_CI = stdp(4) ;
else
    K_i_CI = stdp(1) ;
    K_p_CI = stdp(2) ;
    deltaT_CI = stdp(3) ;
end

% -------------------------------------------------------------------------
%% save to struct
controller_fit_struct.ExprNum = ExprNum ;
controller_fit_struct.MovNum = MovNum ;
controller_fit_struct.K_i = K_i ;
controller_fit_struct.K_p = K_p ;
controller_fit_struct.deltaT = deltaT ;
if KFlag
    controller_fit_struct.K = K ;
end
controller_fit_struct.K_i_CI = K_i_CI ;
controller_fit_struct.K_p_CI = K_p_CI ;
controller_fit_struct.deltaT_CI = deltaT_CI ;
if KFlag
    controller_fit_struct.K_CI = K_CI ;
end
controller_fit_struct.rms = rmse ;
controller_fit_struct.Rsq_adj = Rsq_adj ;
controller_fit_struct.c_pitch = c_pitch ;
%controller_fit_struct.sp_phiR = sp_phiR ;
%controller_fit_struct.sp_phiL = sp_phiL ;
%controller_fit_struct.pitchEstErr = pitchEstErr ;
%controller_fit_struct.phiEstErr = phiEstErr ;
controller_fit_struct.fwdFlipTimes = fwdFlipTimes ;
controller_fit_struct.backFlipTimes = backFlipTimes ;
controller_fit_struct.deltaPhiFront = deltaPhiFront ;
controller_fit_struct.pitchType = pertType ;
controller_fit_struct.rootPath = rootPath ;

if saveFlag
    if ~KFlag
        suffixStr = '_noK_new' ;
    else
        suffixStr = '_LM' ;
    end
    save(fullfile(dataPath, ['controller_fit_struct' suffixStr '.mat']),...
        'controller_fit_struct')
end

end
