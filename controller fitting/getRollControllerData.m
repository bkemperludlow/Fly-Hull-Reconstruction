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
% phir_amp_t = data.phir_amp_t ;
% phil_amp_t = data.phil_amp_t ;
% phiR_amp = data.phiR_amp ;
% phiL_amp = data.phiL_amp ;
bodyRoll = data.anglesLabFrame(:,RHO) ;
t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;
fwdFlipTimesR = data.fwdFlipTimesR ; 
backFlipTimesR = data.backFlipTimesR ; 
fwdFlipTimesL = data.fwdFlipTimesL ; 
backFlipTimesL = data.backFlipTimesL ; 

% match up flip times
[fwdFlipTimesR, fwdFlipTimesL] = alignFlipTimes(fwdFlipTimesR, fwdFlipTimesL) ; 
[backFlipTimesR, backFlipTimesL] = alignFlipTimes(backFlipTimesR, backFlipTimesL) ; 
% ------------------------------------------
%% re-calculate stroke amplitude
% get smoothed wing angles
[~, smoothAnglesR, ~, ~, ~ ] = smoothWingAngles(data, 'R') ;
[~, smoothAnglesL, ~, ~, ~ ] = smoothWingAngles(data, 'L') ;

c_phiR = fit(t', smoothAnglesR(1,:)','cubicinterp') ; 
c_phiL = fit(t', smoothAnglesL(1,:)','cubicinterp') ; 

% -------------------------------------------------------
% calculate stroke amp using smoothed angles (right wing)
phir_t  = [fwdFlipTimesR; backFlipTimesR] ;
phir_atFlip = c_phiR(phir_t) ; 
rmat = [phir_t, phir_atFlip] ; 
rmat = sortrows(rmat,1) ;

phiR_amp = abs(diff(rmat(:,2))) ;
phir_amp_t = rmat(1:end-1,1) + diff(rmat(:,1))/2 ;

% -------------------------------------------------------
% calculate stroke amp using smoothed angles (left wing)
phil_t  = [fwdFlipTimesL; backFlipTimesL] ;
phil_atFlip = c_phiL(phil_t) ; 
lmat = [phil_t, phil_atFlip] ; 
lmat = sortrows(lmat,1) ;

phiL_amp = abs(diff(lmat(:,2))) ;
phil_amp_t = lmat(1:end-1,1) + diff(lmat(:,1))/2 ;

%-------------------------------------------
%% restrict time window
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
% correct for different array lengths
[phir_amp_t, phil_amp_t, R_idx, L_idx] = ...
    alignFlipTimes(phir_amp_t, phil_amp_t) ; 
phiR_amp = phiR_amp(R_idx) ; 
phiL_amp = phiL_amp(L_idx) ;

% average times and get amp difference
phiAmpTimes = (phil_amp_t + phir_amp_t)./2 ; 
phiAmpDiff = phiR_amp - phiL_amp ; 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------
% solve assignment problem for wing flip times
function [timesR_out, timesL_out, R_idx, L_idx] = ...
    alignFlipTimes(timesR, timesL, noMatchCost) 
% if not given, provide cost for not matching (note, here we're using
% seconds)
if ~exist('noMatchCost','var') || isempty(moMatchCost)
    noMatchCost = 1e-3 ; 
end

% make sure the time entries are in column form 
if (size(timesR,1) < size(timesR,2))
    timesR = timesR' ; 
    timesL = timesL' ; 
end

% calculate distances between all cut off times and solve matching problem
distMat = pdist2(timesR, timesL) ; 
M = matchpairs(distMat, noMatchCost) ;
R_idx = M(:,1) ; 
L_idx = M(:,2) ; 
timesR_out = timesR(R_idx) ; 
timesL_out = timesL(L_idx) ; 


end