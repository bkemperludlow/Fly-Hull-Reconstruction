% -------------------------------------------------------------------------
% attempt at a general function for controller fitting -- having so many
% different files for different controller fit types is starting to be a
% headache for making general changes
%
%
% -------------------------------------------------------------------------
function [controller_fit_struct, h] = fitController_dt(rootPath, MovNum, ...
    pertType, pulseTiming, saveFlag, overWriteFlag, plotFlag, debugFlag1, ...
    debugFlag2, smoothFlag, verboseFlag, KFlag)
% --------------------------------
%% inputs
if ~exist('rootPath','var') || isempty(rootPath)
    rootPath = 'D:\Fly Data\Test\14_21102020\'   ; %'D:\Fly Data\VNC Motor Lines\52_19022019\' ;
end
if ~exist('MovNum','var') || isempty(MovNum)
    MovNum = 1 ; %5 ; % 14, 17
end
if ~exist('pertType','var') || isempty(pertType)
    % what type of perturbation (i.e. yaw, pitch, roll) are we fitting?
    pertType = [] ;  % 'Pitch' | 'Yaw' | 'Roll'
end
if ~exist('pulseTiming','var') || isempty(pulseTiming)
    % two element vector giving magnetic pulse start and stop in seconds
    % if not given as input, try to get from experiment info:
    pathStruct = generatePathStruct(rootPath) ; 
    exprInfoStruct = getExprInfo(pathStruct) ; 
    if isfield(exprInfoStruct,'magPulseStart')
        pulseTiming = [exprInfoStruct.magPulseStart, ...
            exprInfoStruct.magPulseEnd] ;
    else
        % if we can't get it, just guess
        pulseTiming = [0, 0.007] ;
    end
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = false  ;
end
if ~exist('overWriteFlag','var') || isempty(overWriteFlag)
    overWriteFlag = false  ;
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true ;
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
if ~exist('KFlag','var') || isempty(KFlag)
    KFlag = false ;
end

% ------------------------------------------------------------------------
%% load fly data
%full data structure
[data, dataPath, pertTypeNum, ExprNum] = ...
    loadPertDataStruct(rootPath, MovNum) ;

% if we can't find data, exit
if isempty(data)
   controller_fit_struct = [] ; 
   h = [] ; 
   return  
end

% get pert type guess from data path
if abs(pertTypeNum) == 1
    pertTypeFromPath = 'Pitch' ;
elseif  abs(pertTypeNum) == 2
    pertTypeFromPath = 'Roll' ;
    % typically need to smooth phiAmpDiff
    smoothFlag = true ;
elseif abs(pertTypeNum) == 3
    pertTypeFromPath = 'Yaw' ;
    % guessing we'll need to smooth deltaAlpha
    smoothFlag = true ;
else
    pertTypeFromPath = [] ; 
end

% see how pert type guess compares to input (fill in if empty)
if isempty(pertType) && ~isempty(pertTypeFromPath)
    pertType = pertTypeFromPath ; 
elseif ~isempty(pertType) && ~isempty(pertTypeFromPath) && ...
        ~strcmp(pertType, pertTypeFromPath)
    fprintf('Warning: mismatch in pert type between input and path \n')
elseif isempty(pertType) && isempty(pertTypeFromPath)
    fprintf('Error determining pert type -- quitting \n')
    controller_fit_struct = [] ;
    return
end

% ------------------------------------------------------
%% set/read fit params
% struct for initial parameter guesses
fit_params = setControllerFitParams(pulseTiming) ;

% read out initial guesses for controller gains
K_i_guess = fit_params.(sprintf('K_i_guess%s',pertType)) ;
K_p_guess = fit_params.(sprintf('K_p_guess%s',pertType)) ;
gainGuess = [K_i_guess, K_p_guess] ; %K_i, K_p, K (from paper)
if KFlag
    K_guess = fit_params.(sprintf('K_guess%s',pertType)) ;
    gainGuess = [gainGuess, K_guess] ;
end

% space of time delays to search
% NB: yaw fits require different time range
if strcmpi(pertType,'Yaw')
    deltaT_min = fit_params.deltaT_minYaw ; %0.002
    deltaT_max =  fit_params.deltaT_maxYaw ; %0.020
    deltaT_N = fit_params.deltaT_NYaw; % 100
else
    deltaT_min = fit_params.deltaT_min ; %0.002
    deltaT_max =  fit_params.deltaT_max ; %0.013
    deltaT_N = fit_params.deltaT_N ; % 100
end

% initial guess for time delay
deltaTGuess = fit_params.(sprintf('deltaTGuess%s',pertType)) ;

% time window in which to fit controller
fitRangeMS = fit_params.(sprintf('%sFitRangeMS',lower(pertType))) ; %  [-10, 45] ; %

% duration/timing of magnetic perturbation
pulseDuration = fit_params.pulseDuration ; % in seconds
pulseTiming = fit_params.pulseTiming ;

% initialize data structure
controller_fit_struct = struct() ;

% ------------------------------------------------------------------------
%% read out relevant data for controller fit
% add some fields to data struct so they can be passed to "get data" func
%data.oneWing = 'R' ;
data.manualCorrRangeMS = fitRangeMS ; % need to pass this to 'getPitchControllerData.m' but it doesn't refer to manual correction any longer
data.pulseDuration = pulseDuration ;
data.pulseTiming = pulseTiming ;

%data to fit controller model to
[wingAngleVals, c_bodyAngle, wingAngleTimes, fwdFlipTimes, ...
    backFlipTimes] = getControllerKinData(data, pertType, debugFlag1,...
    fit_params.steadyPhiFrontFlag) ;

% wingAngleVals = wingAngleVals + 10 ; 

% smooth wing response angles?
if smoothFlag
    wingAngleVals = smooth(wingAngleVals,3)' ;
end

N = length(wingAngleTimes) ;

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

% loop over possible values of time delay
for k = 1:deltaT_N
    % current time delay value
    deltaT = deltaT_space(k) ;
    
    % optimize gain coefficients for current deltaT
    [x, resnorm, residuals, J, Rsq_adj] = ...
        fitControllerGains(wingAngleTimes, wingAngleVals, deltaT, ...
        c_bodyAngle, pertType, gainGuess, false, pulseTiming, KFlag) ;
    
    % store results of current fit
    resNormArray(k) = resnorm ;
    residualArray(k,:) = residuals ;
    K_i_Array(k) = x(1) ;
    K_p_Array(k) = x(2) ;
    Rsq_array(k) = Rsq_adj ;
    J_cell{k} = J ;
    
    % if we solved for constant, store it
    if KFlag
        K_Array(k) = x(3) ;
    end
    
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
    % plot norm of residuals as a function of deltaT?
    figure ;
    hold on
    stem(deltaT_space, resNormArray)
    plot(deltaT_space(bestFitInd), resNormArray(bestFitInd), 'rx')
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
% store common variables shared amongst pert types
controller_fit_struct.ExprNum = ExprNum ; % movie identifiers
controller_fit_struct.MovNum = MovNum ;
controller_fit_struct.pertType = pertTypeNum ;
controller_fit_struct.rootPath = rootPath ;
controller_fit_struct.pulseTiming = pulseTiming ;

controller_fit_struct.K_i = K_i ; % controller params
controller_fit_struct.K_p = K_p ;
controller_fit_struct.deltaT = deltaT ;

controller_fit_struct.K_i_CI = K_i_CI ;  % confidence intervals for controller params
controller_fit_struct.K_p_CI = K_p_CI ;
controller_fit_struct.deltaT_CI = deltaT_CI ;

if KFlag
    controller_fit_struct.K = K ;           % constant term
    controller_fit_struct.K_CI = K_CI ;
end

controller_fit_struct.rms = rmse ;          % fit quality
controller_fit_struct.Rsq_adj = Rsq_adj ;

controller_fit_struct.fwdFlipTimes = fwdFlipTimes ; % wing flip times
controller_fit_struct.backFlipTimes = backFlipTimes ;

% store pert-type specific vars
switch pertType
    case 'Pitch'
        controller_fit_struct.c_pitch = c_bodyAngle ;
        controller_fit_struct.deltaPhiFront = wingAngleVals ;
        
    case 'Roll'
        controller_fit_struct.c_roll = c_bodyAngle ;
        controller_fit_struct.phiAmpTimes = wingAngleTimes ;
        controller_fit_struct.phiAmpDiff = wingAngleVals ;
        
    case 'Yaw'
        controller_fit_struct.c_yaw = c_bodyAngle ;
        controller_fit_struct.midWingBeatTimes = wingAngleTimes ;
        controller_fit_struct.deltaAlpha = wingAngleVals ;
        
    otherwise
        fprintf('Invalid pert type \n')
        keyboard
end

% --------------------------------------------
%% plot final result?
if plotFlag
    [h, ax] = plotControllerFit(controller_fit_struct) ; 
else
    h = [] ; 
end
% --------------------------------------------
%% save results?
if saveFlag
    savePathFull = fullfile(dataPath, 'controller_fit_struct.mat') ;
    
    % check if we should overwrite existing controller fits
    if ~exist(savePathFull, 'file') || overWriteFlag
        save(savePathFull, 'controller_fit_struct')
        
        % also save figure
        if plotFlag
            savePathFig =  fullfile(dataPath, 'controller_fit') ;
            savefig(h, [savePathFig '.fig']) ;
            exportgraphics(h,[savePathFig '.png'], 'Resolution',300) ; 
        end
    else
        % if we don't want to overwrite, spit out a warning. then user can
        % save anyway if they want
        fprintf(['Warning: did not save file: \n %s \n because it' ...
            'already exists \n'], savePathFull)
    end
    
end

end
