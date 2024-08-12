%--------------------------------------------------------------------------
% function to fetch relevant data for PI controller fitting of general DOF
% perturbations
%--------------------------------------------------------------------------
function [wingAngleVals, c_bodyAngle, wingAngleTimes, fwdFlipTimes, ...
    backFlipTimes, steady_t] = getControllerKinData(data, pertType, ...
        debugFlag, steadyPhiFrontFlag)
% -----------------------
%% inputs and params
if ~exist('debugFlag','var')
    debugFlag = false ; 
end
if ~exist('steadyPhiFrontFlag','var')
    steadyPhiFrontFlag = true ; 
end

% general params for getting angle DOF and how to smooth kinematics
defineConstantsScript
smoothingParams = setSmoothingParams ; 

% make angle constant struct to make reading out vals from
% defineConstantsScript easier
angle_constant_struct = struct() ; 
angle_constant_struct.Pitch = THETAB ; 
angle_constant_struct.Roll = RHO ; 
angle_constant_struct.Yaw = PHIB ; 

% other params:
noMatchCost = 1e-3 ; % (sec) used to solve matching problem for flip times
fitSteadyFlag = true ; % if true, find steady-state phi front by linear fit. otherwise find sinle wingstroke with minimal pitch torque
debugFitFlag = false ; % make plots for steady-state estimate?

% pulse timing info
% pulseDuration = data.pulseDuration ; % in milliseconds
pulseTiming = data.pulseTiming ; % in SECONDS ;
pulseTimingMS = 1000*pulseTiming ;  % convert to milliseconds

% determine time range for data that we'll grab
if isfield(data,'manualCorrRangeMS')
    % if we have a range of manually corrected data, or just have already
    % indicated the timing we'd like to use, read that in. however, if
    % manual corr range is too large, we'll need to restrict it
    manualCorrRangeMS_start = max([data.manualCorrRangeMS(1), ...
        (pulseTimingMS(1) - 20)]) ; % -20
    manualCorrRangeMS_end = min([data.manualCorrRangeMS(2), ...
        (pulseTimingMS(1) + 50)]) ;
    manualCorrRangeMS = [manualCorrRangeMS_start, manualCorrRangeMS_end ] ;
else
    % otherwise just take standard guess (although I think in the current
    % version of the code we should never get here)
    manualCorrRangeMS = [pulseTimingMS(1) - 10, pulseTimingMS(1) + 40] ; %[-10, 40] %[-70 70] 
end
manualCorrRange = manualCorrRangeMS / 1000 ; 
%manualCorrRange = [-10,  60] / 1000 ; 

%-----------------------------------------------------
%% read out timing data from struct
fwdFlipTimesR = data.fwdFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
backFlipTimesR = data.backFlipTimesR ;
backFlipTimesL = data.backFlipTimesL ;

t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;

% match up fwd/back flip times
[fwdFlipTimesR, fwdFlipTimesL] = ...
    alignFlipTimes(fwdFlipTimesR, fwdFlipTimesL, noMatchCost) ; 
[backFlipTimesR, backFlipTimesL] = ...
    alignFlipTimes(backFlipTimesR, backFlipTimesL, noMatchCost) ; 

% get left/right mean for forward and back flip times
fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL ) ./ 2 ;
backFlipTimes = (backFlipTimesR + backFlipTimesL) ./ 2 ;

% -------------------------------------------------
%% read out kinematic data (depends on pert type)
% get smoothed wing angles
[~, smoothAnglesR, ~, ~, ~ ] = smoothWingAngles(data, 'R') ;
[~, smoothAnglesL, ~, ~, ~ ] = smoothWingAngles(data, 'L') ;

% get body angle 
bodyAngle = data.anglesLabFrame(:,angle_constant_struct.(pertType)) ;

% filer body angle and get interpolant
filt_lvl = smoothingParams.(sprintf('%s_filt_lvl',lower(pertType))) ;
bodyAngleFilt = filterEulerAngle(bodyAngle, filt_lvl) ; % low-pass butterworth filter

c_bodyAngle = fit(t',bodyAngleFilt,'cubicinterp'); % cubic interpolant

% change which wing angles we grab based on pert type
switch pertType
    case 'Pitch'
        % -------------------
        %% PITCH
        %--------------------
%         % get interpolant for left/right stroke angle 
%         c_phiR = fit(t', smoothAnglesR(1,:)','cubicinterp') ; 
%         c_phiL = fit(t', smoothAnglesL(1,:)','cubicinterp') ; 
%         
        % i know it's a little weird to use spline and interpolant here,
        % but this keeps things matching with older version of code
        phiR = -data.anglesBodyFrame(:,PHIR) ;
        phiL = data.anglesBodyFrame(:,PHIL) ;
        sp_phiR = mySplineSmooth(t(~isnan(phiR)),phiR(~isnan(phiR)),1) ;
        sp_phiL = mySplineSmooth(t(~isnan(phiL)),phiL(~isnan(phiL)),1) ;
        % ---------------------------------------------------------------
        % if we're trying to obtain "phiFrontSteady" (baseline value for
        % phiFront that gives minimal pitch accel), do that:
        if steadyPhiFrontFlag
            [phiFrontSteady, steady_t] = calcSteadyPhiFront(c_bodyAngle,...
                t, fwdFlipTimesR, fwdFlipTimesL, backFlipTimesR, ...
                backFlipTimesL, sp_phiR, sp_phiL,  manualCorrRange, ...
                pulseTiming, fitSteadyFlag, debugFitFlag) ;
        else
            phiFrontSteady = nan ;
            steady_t = nan ;
        end
        
        % ---------------------------------------------------------------
        % get front stroke angles for left and right wing (kind of a pain,
        % but we'll continue to allow the "oneWing" option)
        %% if only using one wing...
        if ~isfield(data,'oneWing')
            phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
            phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;  
        elseif strcmp(data.oneWing,'L')
            phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;  
            phiFrontR = phiFrontL ;
        elseif strcmp(data.oneWing,'R')
            phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
            phiFrontL = phiFrontR ;
        else
            phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
            phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ; 
        end
        
        % ----------------------------------------------
        % get mean phiFront and deltaPhiFront
        phiFront = (phiFrontR + phiFrontL) ./ 2 ;
        
        % use steady phi front value if we have it
        if steadyPhiFrontFlag && ~isnan(phiFrontSteady)
            phiFrontAvg = phiFrontSteady ;
            %disp('worked!')
        else
            prePertFlipIdx = ( fwdFlipTimes <= pulseTiming(1) )  & ...
                ( fwdFlipTimes >= manualCorrRange(1) );
            phiFrontAvg = mean(phiFront(prePertFlipIdx)) ;
        end
        
        % make sure to give final outputs a neutral variable name 
        wingAngleValsR = phiFrontR - phiFrontAvg ; 
        wingAngleValsL = phiFrontL - phiFrontAvg ; 
        wingAngleVals = (phiFront - phiFrontAvg)' ;
        wingAngleTimes = fwdFlipTimes ; 
        
        % label for debug plot
        wingAngleName = '\phi_{fwd}' ; 
        
    case 'Roll'
        % -------------------
        %% ROLL
        %--------------------
        % get interpolant for left/right stroke angle 
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
        % correct for different array lengths
        [phir_amp_t, phil_amp_t, R_idx, L_idx] = ...
            alignFlipTimes(phir_amp_t, phil_amp_t) ;
        phiR_amp = phiR_amp(R_idx) ;
        phiL_amp = phiL_amp(L_idx) ;
        
        % average times and get amp difference
        wingAngleValsR = phiR_amp ; 
        wingAngleValsL = phiL_amp ; 
        wingAngleVals = phiR_amp - phiL_amp ;
        wingAngleTimes = (phil_amp_t + phir_amp_t)./2 ;
        
        % label for debug plot
        wingAngleName = '\Phi_{LR}' ; 
        
    case 'Yaw'
        % -------------------
        %% YAW
        %--------------------
        fprintf(['Warning: this version of the code for yaw' ...
            'perturbations has not been tested! \n'])
        % get wing pitch angles for left/right wing
        etaR = smoothAnglesR(3,:)' ;
        etaL = smoothAnglesL(3,:)' ;
%         
%         c_etaR = fit(t', etaR,'cubicinterp') ; 
%         c_etaL = fit(t', etaL,'cubicinterp') ; 
        
        % -----------------------------------------------------------
        % need to find L/R wing pitch difference at midstroke points
        % From Leif's 2010 paper: "The measured wing orientation angles
        % are used to calculate the attack angle difference Δα, which is 
        % defined to be the right minus left attack angles for the forward 
        % stroke and left minus right angles for the backward stroke."
        % -----------------------------------------------------------
        % combine all wing flip times (both left and right)
        [flipTimesR, flipSortInd] = sort([fwdFlipTimesR ; backFlipTimesR]) ; 
        flipTimesL= sort([fwdFlipTimesL ; backFlipTimesL]) ; 
        
        % get time range in between flip times -- these should correspond
        % to forward and backward stroke
        [~, ~, bin_ind_R] = histcounts(t, flipTimesR) ;
        [~, ~, bin_ind_L] = histcounts(t, flipTimesL) ;
        
        % some time values may fall outside fwd/back stroke bins -- account
        % for these
        zero_bin_idx_R = (bin_ind_R < 1) ;
        zero_bin_idx_L = (bin_ind_L < 1) ;
        
        % get median wing pitch angle during each halfstroke
        alphaR = accumarray(bin_ind_R(~zero_bin_idx_R)', ...
            etaR(~zero_bin_idx_R), [], @nanmedian) ; 
        alphaL = accumarray(bin_ind_L(~zero_bin_idx_L)', ...
            etaL(~zero_bin_idx_L), [], @nanmedian) ; 
        
        % due to Leif's sign convention (see above) need to switch sign for
        % fwd vs back stroke
        if flipSortInd(1) == 1
            % in this case, the first flip time is a forward flip, so the
            % first alpha values correspond to a BACKWARD stroke
            alphaR(2:2:end) = -1.*alphaR(2:2:end) ;
            alphaL(2:2:end) = -1.*alphaL(2:2:end) ;
        else
            % in this case, the first flip time is a back flip, so the
            % first alpha values correspond to a FORWARD stroke
            alphaR(1:2:end) = -1.*alphaR(1:2:end) ;
            alphaL(1:2:end) = -1.*alphaL(1:2:end) ;
        end
        
        % get L/R difference in wing pitch
        wingAngleValsR = alphaR ; 
        wingAngleValsL = alphaL ; 
        wingAngleVals = alphaL - alphaR ; 
        
        midWbTimesR = (flipTimesR(1:end-1) + flipTimesR(2:end))./2 ;
        midWbTimesL = (flipTimesL(1:end-1) + flipTimesL(2:end))./2 ;
        wingAngleTimes = (midWbTimesR + midWbTimesL)./2 ; 
        
        % label for debug plot
        wingAngleName = '\alpha_{LR}' ; 
        
    otherwise
        fprintf('Error: invalid pert type \n')
        keyboard 
end

% ---------------------------------------------------------
%%  restrict time range of wing output to just fit region
% the main thing is to make sure the wing angle var times are restricted:
t_range_idx = (wingAngleTimes >= manualCorrRange(1)) & ...
    (wingAngleTimes <= manualCorrRange(2)) ; 

wingAngleTimes = wingAngleTimes(t_range_idx) ; 
wingAngleVals = wingAngleVals(t_range_idx) ; 
wingAngleValsR = wingAngleValsR(t_range_idx) ; 
wingAngleValsL = wingAngleValsL(t_range_idx) ; 

% since we're outputting fwd and back flip times regardless of pert type,
% restrict range of them too
fwdFlipTimes = fwdFlipTimes((fwdFlipTimes >= manualCorrRange(1)) & ...
    (fwdFlipTimes <= manualCorrRange(2))) ;
backFlipTimes = backFlipTimes((backFlipTimes >= manualCorrRange(1)) & ...
    (backFlipTimes <= manualCorrRange(2))) ;

%--------------------------------------------------------------------------
%% make plot to test data?
if debugFlag
    % -------------------------------------------------
    % get body angle time derivatives over fit range
    t_min = min(wingAngleTimes) ;
    t_max = max(wingAngleTimes) ;
    t_range = linspace(t_min, t_max, 300) ; 
    
    bodyAngleSmoothed = c_bodyAngle(t_range) ;
    [bodyAngleVel, bodyAngleAccel] = differentiate(c_bodyAngle, t_range) ;
    
    % -----------------------------------------
    % make debug figure
    h_debug = figure('Position',[680, 390, 560, 588],...
        'PaperPositionMode','auto') ; 
    
    % ----------------------------
    % body angle subplot
    subplot(3,2,1)
    hold on
    plot(1000*t_range, bodyAngleSmoothed,'b-','LineWidth',1.5)
    plot(1000*t, bodyAngle, 'k.')
    %plot(1000*t_peak, pitch_smoothed(ind_peak), 'ro', 'MarkerFaceColor','r')
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel(sprintf('%s (deg)', pertType))
    xlabel('Time (ms)')
    
    % ----------------------------
    % body vel subplot
    subplot(3,2,3)
    plot(1000*t_range, bodyAngleVel,'r-','LineWidth',1.5)
    ylabel(sprintf('%s Vel (deg/s)', pertType))
    xlabel('Time (ms)')
    set(gca,'xlim', 1000*[t_min t_max])
    
    % ----------------------------
    % body accel subplot
    subplot(3,2,2)
    plot(1000*t_range, bodyAngleAccel,'m-','LineWidth',1.5)
    hold on
    %plot(1000*t_peak, pitchAccel(ind_peak),'bx')
    set(gca,'xlim',1000*[t_min t_max])
    ylabel(sprintf('%s Accel (deg/s^2)', pertType))
    xlabel('Time (ms)')
    
    % -----------------------------
    % wing angle subplot
    subplot(3,2,5)
    hold on
    %plot(1000*t_range, c_phiR(t_range), 'r-')
    %plot(1000*t_range, c_phiL(t_range), 'b-')
    plot(1000*wingAngleTimes, wingAngleVals,'ko-','MarkerFaceColor', 'k')
    %plot(1000*backFlipTimes, phiBack,'k^','MarkerFaceColor', 'k')
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel(sprintf('%s%s (deg)','\Delta', wingAngleName))
    xlabel('Time (ms)')
    
    % -----------------------------
    % wing L/R angle subplot
    subplot(3,2,6)
    hold on
    %plot(1000*t_range, c_phiR(t_range), 'r-')
    %plot(1000*t_range, c_phiL(t_range), 'b-')
    plot(1000*wingAngleTimes, wingAngleValsR,'ro-','MarkerFaceColor', 'r')
    plot(1000*wingAngleTimes, wingAngleValsL,'bo-','MarkerFaceColor', 'b')
    %plot(1000*backFlipTimes, phiBack,'k^','MarkerFaceColor', 'k')
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel(sprintf('L and R %s (deg)', wingAngleName))
    xlabel('Time (ms)')
end

end

