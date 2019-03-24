%--------------------------------------------------------------------------
% function to fetch relevant data for PI controller fitting of roll
% perturbations
%--------------------------------------------------------------------------
function [phiAmpDiff, c_roll, phiAmpTimes] = ...
    getRollControllerData(data, debugFlag)

defineConstantsScript
smoothingParams = setSmoothingParams ;

if isfield(data,'manualCorrRangeMS')
    %manualCorrRangeMS = data.manualCorrRangeMS ;
    manualCorrRangeMS_start = max([data.manualCorrRangeMS(1) -20]) ;
    manualCorrRangeMS_end = min([data.manualCorrRangeMS(2) 60]) ;
    manualCorrRangeMS = [manualCorrRangeMS_start, manualCorrRangeMS_end ] ;
   %manualCorrRangeMS = [-20 55] ;
    %manualCorrRangeMS = [-10 30] ;
else
    manualCorrRangeMS = [-10 50] ;
end
manualCorrRange = manualCorrRangeMS / 1000 ; 

%-------------------------------------------
% read in data
phir_amp_t = data.phir_amp_t ;
phil_amp_t = data.phil_amp_t ;
phiR_amp = data.phiR_amp ;
phiL_amp = data.phiL_amp ;
bodyRoll = data.anglesLabFrame(:,RHO) ;
t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;
%fwdFlipTimesR = data.fwdFlipTimesR ; 
%backFlipTimesR = data.backFlipTimesR ; 

%-------------------------------------------
% restrict time window
correctedIndR = find(phir_amp_t > manualCorrRange(1) & ...
    phir_amp_t < manualCorrRange(2)) ;
correctedIndL = find(phil_amp_t > manualCorrRange(1) & ...
    phil_amp_t < manualCorrRange(2)) ;

phir_amp_t = phir_amp_t(correctedIndR) ;
phil_amp_t = phil_amp_t(correctedIndL) ;
phiR_amp = phiR_amp(correctedIndR) ;
phiL_amp = phiL_amp(correctedIndL) ;

%-------------------------------------------
% low pass filter for body roll angle
roll_filt = filterEulerAngle(bodyRoll,smoothingParams.roll_filt_lvl) ;
c_roll = fit(t',roll_filt,'cubicinterp');

%-------------------------------------------
% this doesn't really make sense here...
if isfield(data,'oneWing') && strcmp(data.oneWing,'L')
    phiR_amp = phiL_amp ; 
    phir_amp_t = phil_amp_t ;
elseif isfield(data,'oneWing') && strcmp(data.oneWing,'R')
    phiL_amp = phiR_amp ; 
    phil_amp_t = phir_amp_t ;
end
    
%-------------------------------------------
% correct for different array lengths
if length(phir_amp_t) == length(phil_amp_t)
    phiAmpDiff = phiR_amp - phiL_amp ;
    phiAmpTimes = (phil_amp_t + phir_amp_t ) /2 ;
elseif length(phir_amp_t) < length(phil_amp_t)
    idx = zeros(length(phir_amp_t),1) ;
    for q = 1:length(phir_amp_t)
        [~,minInd] = min(abs(phil_amp_t - phir_amp_t(q))) ; 
        idx(q) = minInd ;
    end
    phiAmpTimes = (phir_amp_t + phil_amp_t(idx)) / 2 ;
    phiAmpDiff = (phiR_amp - phiL_amp(idx)) / 2 ;
elseif length(phil_amp_t) < length(phir_amp_t)
    idx = zeros(length(phil_amp_t),1) ;
    for q = 1:length(phil_amp_t)
        [~,minInd] = min(abs(phir_amp_t - phil_amp_t(q))) ; 
        idx(q) = minInd ;
    end
    phiAmpTimes = (phir_amp_t(idx) + phil_amp_t) / 2 ;
    phiAmpDiff = (phiR_amp(idx) - phiL_amp) / 2 ;
end

%--------------------------------------------------------------------------
% make plot to test data?
if debugFlag
    xlim = 1000*[manualCorrRange(1), max(phiAmpTimes)] ; 
    h_debug = figure('Position',[680, 390, 560, 588],'PaperPositionMode','auto') ; 
    
    subplot(3,1,1)
    hold on
    plot(phir_amp_t*1000, phiR_amp,'ro-')
    plot(phil_amp_t*1000, phiL_amp,'bo-')
    set(gca,'xlim',xlim)
    ylabel({'Stroke Amp.', '\Phi , [deg]'})
    legend({'\Phi_R','\Phi_L'},'location', 'northwest')
    grid on
    
    subplot(3,1,2)
    plot(phiAmpTimes*1000, phiAmpDiff ,'ko-', 'markerfacecolor', 'k')
    %xlabel('Time [ms]')
    ylabel({'Stroke Amp. Diff.', '\Delta\Phi , [deg]'})
    set(gca,'xlim',xlim)
    grid on

    subplot(3,1,3)
    hold on
    plot(t*1000, bodyRoll, 'k.')
    plot(t*1000, roll_filt, 'r-')
    %xlabel('Time [ms]')
    ylabel({'Body Roll', '\rho , [deg]'})
    set(gca,'xlim',xlim)
    %set(gca,'xlim',[phiAmpTimes(1)*1000  phiAmpTimes(end)*1000])
    max_rho = max(bodyRoll((t > phiAmpTimes(1) & t < phiAmpTimes(end)))) ; 
    min_rho = min(bodyRoll((t > phiAmpTimes(1) & t < phiAmpTimes(end)))) ; 
    set(gca,'ylim',[(min_rho-5) (max_rho+5)]) ; 
    grid on
    %keyboard ;
end

phiAmpTimes = phiAmpTimes' ; 
phiAmpDiff = phiAmpDiff' ; 

end