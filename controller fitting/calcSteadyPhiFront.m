% -------------------------------------------------------------------------
% function to estimate the baseline forward stroke amplitude that leads to
% zero pitching torque -- this is used as the reference value for pitch
% controller fits.
%
% attempts to do this by performing a linear fit to pitch acceleration vs
% phi front, and locates the value of pitch_accel = 0. Failing that, tries
% to find the value for phi front that has minimal absolute pitch accel
%
% -------------------------------------------------------------------------
function [phiFrontSteady, steady_t] = calcSteadyPhiFront(c_pitch, t,...
    fwdFlipTimesR, fwdFlipTimesL, backFlipTimesR, backFlipTimesL, ...
    c_phiR, c_phiL,  manualCorrRange, pulseTiming, fitSteadyFlag, debugFlag)
% -------------------------------
% inputs
if ~exist('fitSteadyFlag','var') || isempty(fitSteadyFlag)
   fitSteadyFlag = true ;  
end
if ~exist('debugFlag','var') || isempty(debugFlag)
   debugFlag = false ;  
end

% ----------------------------------------------
% params
fitThreshRMS_high = 8 ; % max value for deltaPhiFront RMS from fit in general
fitThreshRMS_low = 3 ; % max value for deltaPhiFront RMS from fit if we have few points to fit
fitMinDOF = 15 ;   % minimum number of points we need to take fit as good
fitRsqThresh = 0 ; % min rsq value. currently set to zero, so i guess it's not doing anything atm

% --------------------------------------------
% calculate body pitch acceleration
[~, pitch_accel] = differentiate(c_pitch,t) ;

% also get fwd flip positions (left, right, and mean)
% phiFrontTempR = c_phiR(fwdFlipTimesR) ;
% phiFrontTempL = c_phiL(fwdFlipTimesL) ;
phiFrontTempR = fnval(c_phiR,fwdFlipTimesR) ;
phiFrontTempL = fnval(c_phiL,fwdFlipTimesL) ;
phiFrontTempMean = (phiFrontTempR + phiFrontTempL)./2 ;

% --------------------------------------------------------------------
% want to get average pitch acceleration for each wingstroke 
% (binned between back flip times). so first get mean back flip times
meanBackFlipTimes = (backFlipTimesR + backFlipTimesL) ./ 2 ;
binCenters = (meanBackFlipTimes(2:end) + meanBackFlipTimes(1:end-1))./2 ;
% N_bins = length(meanBackFlipTimes)-1 ;

% loop through bins
[~, ~, bin_ind] = histcounts(t, meanBackFlipTimes) ;
zero_bin_idx = (bin_ind < 1) ; 
binnedPitchAccel = accumarray(bin_ind(~zero_bin_idx)', ...
    pitch_accel(~zero_bin_idx), [], @nanmean) ; 

% ------------------------------------------------------------------
% find wingbeat with least pitch acceleration
if fitSteadyFlag
    % this version tries to fit relationship between phi front and
    % pitch  accel (steady value is when linear fit = 0)
    
    % find pert time to exclude from fit
    pert_idx = (binCenters >= pulseTiming(1)) & ...
        (binCenters <= pulseTiming(2)) ;
    % make sure fwd flip times line up with the right bin center for the
    % acceleration (nb: bin centers are halfway point between back flip
    % points, so they should be close to fwd flip points)
    [~, ~, fwdInd, ~] = alignFlipTimes(fwdFlipTimesR, binCenters) ;
    
    % do linear fit for pitch accel vs phi front
    if length(fwdInd) == length(pert_idx)
        [c_phi_accel, gof, output] = fit(phiFrontTempMean(fwdInd & ~pert_idx), ...  % phiFrontTempMean(fwdInd & ~pert_idx)
            binnedPitchAccel(~pert_idx),'poly1','Robust','LAR') ;
    else
        % getting some weird error here; should try to fix later
        phiFrontSteady = nan ;
        steady_t = nan ; 
        return
    end
    % take phiFrontSteady as the value of phi front that gives 0 pitch
    % accel (according to linear fit)
    phiFrontSteady = (-1*c_phi_accel.p2 / c_phi_accel.p1) ;
    
    % ---------------------------------------------------------------------
    % make sure this is a reasonable value -- look at what this would give
    % for deltaPhiFront in time prior to pert
    % ---------------------------------------------------------------------
    % pre-perturbation index
    pre_pert_idx = (fwdFlipTimesR >= manualCorrRange(1)) & ...
        (fwdFlipTimesR <= pulseTiming(1));
    
    % delta phi front (according to set point from fit) in pre-pert time
    % range and its RMS value
    testDeltaPhiFront = phiFrontTempMean(pre_pert_idx) - phiFrontSteady ;
    testDeltaPhiFrontRMS = sqrt(mean((testDeltaPhiFront).^2)) ;
    
    % if we're fitting to few points, require very low rms deltaPhiFront
    dof_check = (testDeltaPhiFrontRMS >= fitThreshRMS_low) & ...
        (gof.dfe < fitMinDOF) ; 
    % in general, don't want to have high rms deltaPhiFront
    rms_check = testDeltaPhiFrontRMS >= fitThreshRMS_high ; 
    % want rsquare to be greater than zero
    rsq_check = gof.rsquare < fitRsqThresh ; 
    % make sure solver ran its way through
    exitflag_check = output.exitflag < 1 ; 
    % make sure slope of fit is negative (based on our sign convention)
    slope_check = c_phi_accel.p1 > 0 ;
    
    % if any of the above checks are true, just return nan (then we'll
    % estimate deltaPhiFront by subtracting mean of pre-pert values)
    if dof_check || rms_check || rsq_check || exitflag_check || slope_check
        phiFrontSteady = nan ;
    end
    
    % not sure why i had this in previous version, but for backwards
    % compatibility, need to assign a value to "steady_t" (which is the
    % time value of the front stroke with minimum absolute pitch accel).
    % since this doesn't make sense in the context of the fit strategy,
    % just return nan
    steady_t = nan ; 
    
else
    % this version finds the wingbeat (prior to pert) with lowest
    % pitch acceleration
    [~, sortIdx] = sort(abs(binnedPitchAccel),'ascend') ;
    binTimesSorted = binCenters(sortIdx) ;
    steadyBinNum = find(binTimesSorted < 0, 1, 'first') ;
    steadyBinTime = binTimesSorted(steadyBinNum) ;
    
    % find the forward stroke angle that this corresponds to
    [~, steadyTimeInd] = min(abs(t - steadyBinTime)) ;
    steady_t = t(steadyTimeInd) ;
    [~, steadyIndR] = min(abs(fwdFlipTimesR - steady_t)) ;
    [~, steadyIndL] = min(abs(fwdFlipTimesL - steady_t)) ;
    phiFrontSteady = (phiFrontTempR(steadyIndR) + ...
        phiFrontTempL(steadyIndL))./2  ;
end

% ----------------------------------------------------------------
% debug?
if debugFlag
    figure ;
    hold on
    if fitSteadyFlag
        plot(phiFrontTempR(fwdInd), binnedPitchAccel,'ro')
        plot(phiFrontTempL(fwdInd), binnedPitchAccel,'bo')
        plot(phiFrontTempMean(fwdInd), binnedPitchAccel,'ko')
        x_temp = linspace(min(phiFrontTempMean), max(phiFrontTempMean), 20) ;
        plot(x_temp, c_phi_accel(x_temp),'k-')
        xlabel('\Delta \phi_{front} (deg)')
        ylabel('Pitch Accel (deg/s^2)')
        title('Steady state \phi_{front} fit')
    else
        plot(t, pitch_accel)
        axis tight
        plot(binCenters, binnedPitchAccel,'ko-')
        plot(steadyBinTime, binnedPitchAccel(sortIdx(steadyBinNum)),'rx')
        
        xlabel('Time (s)')
        ylabel('Pitch Accel (deg/s^2)')
        title('Wingbeat-averaged pitch accel')
    end
    keyboard
end

end