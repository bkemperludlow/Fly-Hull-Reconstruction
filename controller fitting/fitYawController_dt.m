

rootPath = 'I:\B1 Yaw Data\Yaw analysis\Saumya Manually Corrected\18_02062017\Analysis\' ;
%rootPath = 'H:\Fly Data\Janelia Flies\kir round 2\19_31102016\Analysis\' ;
ExprNum = 18 ;
MovNum = 58 ;
pertType = nan ;
flyType = 2 ;
plotFlag1 = false ;
plotFlag2 = true ;
debugFlag1 = false ;
debugFlag2 = true ;
saveFlag = false ; 

gainGuess = [0, 0.008] ; %K_i, K_p, K (from paper)

% controller_fit_struct = struct('ExprNum',[],'MovNum',[],'K_i',[],'K_p',[],...
%     'deltaT',[],'K',[],'K_i_CI',[],'K_p_CI', [], 'K_CI', [], 'deltaT_CI', [],...
%     'rms',[],'c_pitch',[],'phiEstErr',[],'deltaPhiFront',[],'fwdFlipTimes',[],...
%     'backFlipTimes', [],'t',[],'flyType',[]) ;
controller_fit_struct = struct() ; 

%% get data

%full data structure
data = loadPertDataStruct(rootPath, ExprNum, MovNum, pertType) ;
%data.oneWing = 'L' ; 
%data.manualCorrRangeMS = [-10 , 40] ;

%data to fit controller model to
[deltaAlphaMean, c_yaw, midWingBeatTimes, deltaAlphaSEM] = ...
                                   getYawControllerData(data, debugFlag1) ;

N = length(midWingBeatTimes) ;

%% make discretized t space

deltaT_min = 0.007 ;
deltaT_max = 0.025 ;
deltaT_N = 80 ;

deltaT_space = linspace(deltaT_min, deltaT_max, deltaT_N) ;

%% loop through deltaT and fit gains
resNormArray = nan(deltaT_N,1) ;
residualArray = nan(deltaT_N,N) ;
K_i_Array = nan(deltaT_N,1) ;
K_p_Array = nan(deltaT_N,1) ;
J_cell = cell(deltaT_N,1) ;

%cc = 1 ;
for k = 1:deltaT_N
    deltaT = deltaT_space(k) ;
    [K_i, K_p, resnorm, residuals, J] = ...
        fitYawControllerGains(midWingBeatTimes, deltaAlphaMean, deltaT,...
        c_yaw, gainGuess, plotFlag1)  ;
    if plotFlag1
        title(['\DeltaT = ' num2str(deltaT)])
    end
    
    resNormArray(k) = resnorm ;
    residualArray(k,:) = residuals ;
    K_i_Array(k) = K_i ;
    K_p_Array(k) = K_p ;
    J_cell{k} = J ;
    
    %cc = cc + 1 ;
end

%% estimate jacobian for deltaT
J_deltaT_Array = zeros(N, deltaT_N) ;

for i = 1:N
    residuals = residualArray(:,i) ;
    c_residual = fit(deltaT_space', residuals, 'cubicinterp') ;
    
    residual_derivative = differentiate(c_residual, deltaT_space) ;
    J_deltaT_Array(i,:) = residual_derivative' ;
end
%% find best fit

[~, bestFitInd] = min(resNormArray) ;
K_i = K_i_Array(bestFitInd) ;
K_p = K_p_Array(bestFitInd) ;
deltaT = deltaT_space(bestFitInd) ;
resnorm = resNormArray(bestFitInd) ;
rmse = sqrt(resnorm/N) ;

Jacobian_gains = J_cell{bestFitInd} ;
Jacobian_gains = full(Jacobian_gains) ;

J_deltaT = J_deltaT_Array(:, bestFitInd) ;
Jacobian = [Jacobian_gains, J_deltaT] ;

varp = resnorm*inv(Jacobian'*Jacobian)/N ;
stdp = sqrt(diag(varp));
%stdp = 100*stdp'./p;

if debugFlag2
    figure ;
    hold on
    stem(deltaT_space, resNormArray)
    plot(deltaT_space(bestFitInd), resNormArray(bestFitInd), 'rx')
end

if plotFlag2
    [~,~, ~, ~, ~] = ...
        fitYawControllerGains(midWingBeatTimes, deltaAlphaMean, ...
            deltaT_space(bestFitInd), c_yaw, gainGuess, true)  ;
end

disp(['K_i = ' num2str(K_i) ' +/- ' num2str(stdp(1))])
disp(['K_p = ' num2str(K_p) ' +/- ' num2str(stdp(2))])
disp(['DT = ' num2str(deltaT) ' +/- ' num2str(stdp(3))])

disp('RMSE:')
disp(rmse)

K_i_CI = stdp(1) ; K_p_CI = stdp(2) ; deltaT_CI = stdp(3) ; 

if flyType == 1
    plotColor = [.7 0 0 ] ;
    flyTypeStr = 'experimental' ;
elseif flyType == 2
    plotColor = [0 .7 0 ] ;
    flyTypeStr = 'control' ;
else
    plotColor = .5*[1 1 1] ;
end

%% save to struct
controller_fit_struct.ExprNum = ExprNum ;
controller_fit_struct.MovNum = MovNum ;
controller_fit_struct.K_i = K_i ;
controller_fit_struct.K_p = K_p ;
controller_fit_struct.deltaT = deltaT ;
controller_fit_struct.K = [] ;
controller_fit_struct.K_i_CI = K_i_CI ;
controller_fit_struct.K_p_CI = K_p_CI ;
controller_fit_struct.deltaT_CI = deltaT_CI ;
controller_fit_struct.K_CI =[] ;
controller_fit_struct.rms = rmse ;
controller_fit_struct.flyType = flyTypeStr ;
controller_fit_struct.c_pitch = c_yaw ;
%controller_fit_struct.sp_phiR = sp_phiR ;
%controller_fit_struct.sp_phiL = sp_phiL ;
%controller_fit_struct.pitchEstErr = pitchEstErr ;
%controller_fit_struct.phiEstErr = phiEstErr ;
controller_fit_struct.midWingBeatTimes = midWingBeatTimes ;
controller_fit_struct.deltaAlphaMean = deltaAlphaMean ;
controller_fit_struct.deltaAlphaSEM = deltaAlphaSEM ;
controller_fit_struct.pitchType = pertType ;

if saveFlag
    save controller_fit_struct_LM controller_fit_struct
end

