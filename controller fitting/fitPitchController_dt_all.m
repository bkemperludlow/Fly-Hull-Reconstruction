logPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis' ;
currPath = pwd ; 
cd(logPath)
rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ;
%% set up loop 

dataLog = importdata('Pitch_Controller_Log.xlsx') ;
ExprNumArray = dataLog.data.Finished(:,1) ;
MovNumArray = dataLog.data.Finished(:,2) ;
flyTypeArray = dataLog.data.Finished(:,5) ;
pitchTypeArray = dataLog.data.Finished(:,3) ;

ExprNumArray = ExprNumArray(~isnan(ExprNumArray)) ;
MovNumArray = MovNumArray(~isnan(MovNumArray)) ;
flyTypeArray = flyTypeArray(~isnan(flyTypeArray)) ;
pitchTypeArray = pitchTypeArray(~isnan(pitchTypeArray)) ;

%rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ;
%ExprNum = 1 ;
%MovNum = 52 ;
%pitchType = -1 ;
plotFlag1 = false ;
plotFlag2 = false ;
debugFlag1 = false ;
debugFlag2 = false ;

saveFlag = true ; 

gainGuess = [0, 0.008, 0] ; %K_i, K_p, K (from paper)
gainGuess_P = [0.008, 0] ; 
cc = 1 ; 

% controller_fit_struct = struct('ExprNum',[],'MovNum',[],'K_i',[],'K_p',[],...
%         'deltaT',[],'K',[],'K_i_CI',[],'K_p_CI', [], 'K_CI', [], 'deltaT_CI', [],...
%         'rms',[],'c_pitch',[],'phiEstErr',[],'deltaPhiFront',[],'fwdFlipTimes',[],...
%         'backFlipTimes', [],'t',[],'flyType',[]) ;

controller_fit_struct = struct() ; 
    
for i = 1:length(MovNumArray)
    
    
%% get data
    ExprNum = ExprNumArray(i) ;
    MovNum = MovNumArray(i) ;
    flyType = flyTypeArray(i) ;
    pitchType = pitchTypeArray(i) ; 
    %full data structure
    data = loadPertDataStruct(rootPath, ExprNum, MovNum, pitchType) ;
    
    %data to fit controller model to
    [deltaPhiFront, c_pitch,fwdFlipTimes, backFlipTimes] = getPitchControllerData(data, debugFlag1) ;
    
    N = length(fwdFlipTimes) ;
    %N = length(fwdFlipTimes(1): 0.125e-3 : fwdFlipTimes(end)) ;
    %% make discretized t space
    
    deltaT_min = 0.004 ;
    if ExprNum == 11 && MovNum == 33
        deltaT_max = 0.01 ;
    else
        deltaT_max = 0.017 ;
    end
    deltaT_N = 60 ;
    
    deltaT_space = linspace(deltaT_min, deltaT_max, deltaT_N) ;
    
    %% loop through deltaT and fit gains
    resNormArray = nan(deltaT_N,1) ;
    residualArray = nan(deltaT_N,N) ;
    K_i_Array = nan(deltaT_N,1) ;
    K_p_Array = nan(deltaT_N,1) ;
    K_Array = nan(deltaT_N,1) ;
    J_cell = cell(deltaT_N,1) ;
    
    %cc = 1 ;
    for k = 1:deltaT_N
        deltaT = deltaT_space(k) ;
        [K_i, K_p, K, resnorm, residuals, J] = ...
            fitPitchControllerGains(fwdFlipTimes, deltaPhiFront, deltaT, c_pitch, gainGuess, plotFlag1)  ;
        if plotFlag1
            title(['\DeltaT = ' num2str(deltaT)])
        end
        
        resNormArray(k) = resnorm ;
        residualArray(k,:) = residuals ;
        K_i_Array(k) = K_i ;
        K_p_Array(k) = K_p ;
        K_Array(k) = K ;
        J_cell{k} = J ;
        
        %cc = cc + 1 ;
    end
    
    %% estimate jacobian for deltaT
    J_deltaT_Array = zeros(N, deltaT_N) ;
    
    for j = 1:N
        residuals = residualArray(:,j) ;
        c_residual = fit(deltaT_space', residuals, 'cubicinterp') ;
        
        residual_derivative = differentiate(c_residual, deltaT_space) ;
        J_deltaT_Array(j,:) = residual_derivative' ;
    end
    %% find best fit
    
    [~, bestFitInd] = min(resNormArray) ;
    K_i = K_i_Array(bestFitInd) ;
    K_p = K_p_Array(bestFitInd) ;
    K = K_Array(bestFitInd) ;
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
         xlabel('\Delta T')
         ylabel('RMS')
         
        figure ; 
         scatter(deltaT_space,K_i_Array,resNormArray/(0.01*max(resNormArray)),...
             resNormArray,'filled')
         xlabel('\Delta T')
         ylabel('K_i')
    end
    
    if plotFlag2
        [~,~, ~, ~, ~, ~] = ...
            fitPitchControllerGains(fwdFlipTimes, deltaPhiFront, deltaT_space(bestFitInd), c_pitch, gainGuess, true)  ;
    end
    
    disp(['K_i = ' num2str(K_i) ' +/- ' num2str(stdp(1))])
    disp(['K_p = ' num2str(K_p) ' +/- ' num2str(stdp(2))])
    disp(['K = ' num2str(K) ' +/- ' num2str(stdp(3))])
    disp(['DT = ' num2str(deltaT) ' +/- ' num2str(stdp(4))])
    
    disp('RMSE:')
    disp(rmse)
    
    K_i_CI = stdp(1) ; K_p_CI = stdp(2) ; K_CI = stdp(3) ; deltaT_CI = stdp(4) ;
    
    if flyType == 1
        plotColor = [.7 0 0 ] ;
        flyTypeStr = 'experimental' ;
    elseif flyType == 2
        plotColor = [0 .7 0 ] ;
        flyTypeStr = 'control' ;
    else
        plotColor = .5*[1 1 1] ;
    end

    controller_fit_struct(cc).ExprNum = ExprNum ;
    controller_fit_struct(cc).MovNum = MovNum ;
    controller_fit_struct(cc).K_i = K_i ;
    controller_fit_struct(cc).K_p = K_p ;
    controller_fit_struct(cc).deltaT = deltaT ;
    controller_fit_struct(cc).K = K ;
    controller_fit_struct(cc).K_i_CI = K_i_CI ;
    controller_fit_struct(cc).K_p_CI = K_p_CI ;
    controller_fit_struct(cc).deltaT_CI = deltaT_CI ;
    controller_fit_struct(cc).K_CI = K_CI ;
    controller_fit_struct(cc).rms = rmse ;
    controller_fit_struct(cc).flyType = flyTypeStr ;
    controller_fit_struct(cc).c_pitch = c_pitch ;
    %controller_fit_struct.sp_phiR = sp_phiR ;
    %controller_fit_struct.sp_phiL = sp_phiL ;
    %controller_fit_struct.pitchEstErr = pitchEstErr ;
    %controller_fit_struct.phiEstErr = phiEstErr ;
    controller_fit_struct(cc).fwdFlipTimes = fwdFlipTimes ;
    controller_fit_struct(cc).backFlipTimes = backFlipTimes ;
    controller_fit_struct(cc).deltaPhiFront = deltaPhiFront ;
    controller_fit_struct(cc).pitchType = pitchType ;
    
    %save controller_fit_struct_LM controller_fit_struct
    cc = cc + 1 ; 
end

if saveFlag 
    save(fullfile(logPath,'pitch_controller_fits_all.mat'),'controller_fit_struct')
end

cd(currPath)


%{


exp_ind = arrayfun(@(x) strcmp(x.flyType,'experimental'),controller_fit_struct) ;
ctrl_ind = arrayfun(@(x) strcmp(x.flyType,'control'),controller_fit_struct) ;
t_bound_ind = ([controller_fit_struct.deltaT] < 0.00401) ; 

Ki_exp = [controller_fit_struct(exp_ind & ~t_bound_ind).K_i] ;
Ki_ctrl = [controller_fit_struct(ctrl_ind & ~t_bound_ind).K_i] ;

Kp_exp = [controller_fit_struct(exp_ind).K_p] ;
Kp_ctrl = [controller_fit_struct(ctrl_ind).K_p] ;

figure ; hold on
plot(ones(size(Ki_exp)),Ki_exp,'rx')
plot(2*ones(size(Ki_ctrl)),Ki_ctrl,'bo')

set(gca,'xlim',[0.5,2.5])

%}