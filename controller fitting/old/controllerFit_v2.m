function [results_full, CI, PitchEstErr] = ...
    controllerFit_v2(ExprNum,MovNum,PitchType,varyPitchSplineFlag,...
        plotFlag,plotFlag2,plotFlag3)
%{
function [results_full, results_NoK, results_NoKi, results_SameT,...
    results_SameT_NoK, results_NoKp, results_NoKp_NoK, PitchEstErr] = ...
    controllerFit_v2(ExprNum,MovNum,PitchType,varyPitchSplineFlag,...
        plotFlag,plotFlag2,plotFlag3)
%}
%Create 3D space of control parameters, find minimum

defineConstantsScript
patchColor = [1 1 1]*.8 ;
%{
ExprNum = 7 ; 
MovNum = 9 ; 
PitchType = 'up' ;
plotFlag = false ; %check the splines
varyPitchSplineFlag = false ;
plotFlag2 = true ; %visualize controller phase space with contour plots
plotFlag3 = true ; %Show control fit with data
%}
%% Load data

if MovNum < 10
    zstr = '00';
elseif MovNum < 100
    zstr = '0';
else
    zstr = '';
end

if strcmp(PitchType,'up')
    datapath = strcat('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Up\Expr_',...
        num2str(ExprNum), '_mov_',zstr,num2str(MovNum)) ;
elseif strcmp(PitchType,'down')
    datapath = strcat('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Down\Expr_',...
        num2str(ExprNum), '_mov_',zstr,num2str(MovNum)) ;
else
    disp('check PitchType')
    return ;
end
cd(datapath)

datafilename = strcat(datapath,'\Expr',num2str(ExprNum), ...
    'mov',zstr,num2str(MovNum), '_Data_manually_corrected_onlyPhiFront.mat') ;

load(datafilename) ;

%Check ignoreFrames and correctionTime
if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

if (isfield(data,'correctionTime'))
    correctionTime=data.correctionTime;
elseif(isfield(data,'manualCorrRangeMS'))
    correctionTime=data.manualCorrRangeMS;
end

if isfield(data.params,'pulseLengthMS')
    pulseLength = data.params.pulseLengthMS ;
elseif (~isfield(data.params,'pulseLengthMS')) && (ExprNum == 7)
    pulseLength = 5.8 ; %ms
else
    pulseLength = 7 ; %ms
end

%Get flip times
if isfield(data,'fwdFlipTimesR') && isfield(data,'backFlipTimesR')
    fwdFlipTimesR = data.fwdFlipTimesR ;
    backFlipTimesR = data.backFlipTimesR ;
    fwdFlipTimesL = data.fwdFlipTimesL ;
    %backFlipTimesL = data.backFlipTimesL ;
else
    [fwdFlipTimesR, ~, fwdFlipTimesL, ~, ~,~, data]...
        = saveWingFlipsAndAngles(ExprNum,MovNum,PitchType) ;
end

%Get the body pitch angle
if (isfield(data, 'anglesLabFrame'))
    bodyPitch = data.anglesLabFrame(:,BETA) ;
else
    [~, ~, ~, ~, ~, anglesLabFrame, data] = ...
        saveWingFlipsAndAngles(ExprNum,MovNum,PitchType) ;
    bodyPitch = anglesLabFrame(:,BETA) ;
end

%Get wing stroke angle
phiR = (180/pi)*unwrap(-(pi/180)*data.anglesBodyFrame(:,PHIR)) ; 
phiL = (180/pi)*unwrap((pi/180)*data.anglesBodyFrame(:,PHIL)) ; 

%% Exclude bad frames and smooth
startTime = data.params.startTrackingTime ;
endTime = data.params.endTrackingTime ;
fps = data.params.fps ;
t = (startTime:endTime)/fps ; %seconds
tms = 1000*t ; %ms
dt = 1/fps ;
[~,correctioStartInd] = min(abs(t - correctionTime(1)/1000 + .025)) ;%+.025)) ;
zeropoint = find(t == 0) ;

bodyPitch(ignoreFrames) = NaN ;

%Find bad phi frames (there are more)
badPhiRFrames = find((phiR > 185) | (phiR < 0)) ; 
badPhiLFrames = find((phiL > 185) | (phiL < 0)) ; 

phiRInd = [badPhiRFrames',ignoreFrames];
phiR(phiRInd) = NaN ;
phiLInd = [badPhiLFrames',ignoreFrames];
phiL(phiLInd) = NaN ; 

%Now look at only frames that are manually corrected (AC = 'after correction')
t_AC = t(correctioStartInd:end) ; 
bodyPitch_AC = bodyPitch(correctioStartInd:end) ;
phiR_AC = phiR(correctioStartInd:end) ;
phiL_AC = phiL(correctioStartInd:end) ; 

%Smooth pitch
ind = find(~isnan(bodyPitch_AC)) ;
currtvec = t_AC(ind) ;
currBodyPitch = bodyPitch_AC(ind)' ;

if varyPitchSplineFlag
    PitchErrLow = .3 ;
    PitchErrHigh = .8 ;
    Ntest = 30 ;
    maxRMS = 1 ;
    
    [sp_pitch, pitch_smooth, RMS, PitchEstErr] = varySpline(currtvec,currBodyPitch,...
        PitchErrLow,PitchErrHigh,Ntest,maxRMS,plotFlag); 
else
    PitchEstErr = .6 ;
    [sp_pitch, pitch_smooth, ~] = mySplineSmooth(currtvec, currBodyPitch, PitchEstErr) ;
end

%Calculate velocity spline and find average pitch to subtract off
sp_velocity = fnder(sp_pitch,1); 
bodyPitchAvg = nanmean(bodyPitch(correctioStartInd:zeropoint)) ;

%Smooth phiR
ind2 = find(~isnan(phiR_AC)) ; 
currPhiR = phiR_AC(ind2) ; 
currtvec2 = t_AC(ind2) ; 

PhiREstErr = 1.5 ;
[sp_phiR, phiR_smooth, ~] = mySplineSmooth(currtvec2,currPhiR,PhiREstErr) ;

%Smooth phiL
ind3 = find(~isnan(phiL_AC)) ; 
currPhiL = phiL_AC(ind3) ; 
currtvec3 = t_AC(ind3) ; 

PhiLEstErr = 1.5 ;
[sp_phiL, phiL_smooth, ~] = mySplineSmooth(currtvec3,currPhiL,PhiLEstErr) ;

%% Check with plots (if necessary)
if plotFlag
   figure ; %pitch angle
   plot(t_AC*1000, fnval(sp_pitch,t_AC), 'g-')
   hold on 
   plot(t_AC*1000, bodyPitch_AC,'k.')
   
   figure ; %pitch velocity
   plot(t_AC*1000, fnval( fnder(sp_pitch,1),t_AC), 'k-')
   
   figure ; %right wingstroke angle
   plot(t_AC*1000, fnval(sp_phiR,t_AC), 'r-')
   hold on 
   plot(t_AC*1000, phiR_AC,'ro')
   
   figure ; %left wingstroke angle
   plot(t_AC*1000, fnval(sp_phiL,t_AC), 'b-')
   hold on 
   plot(t_AC*1000, phiL_AC,'bo')
end

%% Initialize control space

ki_min = -1 ; %unitless
ki_max = 1 ;
ki_N = 25 ;

kp_min = 0 ; %seconds
kp_max = .05 ; 
kp_N = 25 ;

deltat_min = .004 ; %seconds
deltat_max = .03 ;
deltat_N = 25 ; 

%deltat2_min = .004 ; %seconds
%deltat2_max = .03 ;
%deltat2_N = 30 ;

k_min = -4 ;
k_max = 4 ;
k_N = 9 ; 

ki_lin = linspace(ki_min, ki_max, ki_N) ;
kp_lin = linspace(kp_min, kp_max, kp_N) ;
deltat_lin = linspace(deltat_min, deltat_max, deltat_N) ;
%deltat2_lin = linspace(deltat2_min, deltat2_max, deltat2_N) ;
k_lin = linspace(k_min, k_max, k_N) ; 

systemSize = ki_N*kp_N*deltat_N*k_N ;

[Ki,Kp,DeltaT,K] = ndgrid(ki_lin,kp_lin,deltat_lin,k_lin) ; 
dim = size(Ki) ;
ErrorList = nan(systemSize,1) ; 

%% Find values of phi_front after perturbation (or correction)


%firstFwdIndR = find(fwdFlipTimesR > pulseLength/1000, 1, 'first') ;
%firstFwdIndL = find(fwdFlipTimesL > pulseLength/1000, 1, 'first') ;

firstFwdIndR = find(fwdFlipTimesR > correctionTime(1)/1000, 1, 'first') ;
firstFwdIndL = find(fwdFlipTimesL > correctionTime(1)/1000, 1, 'first') ;

%Check to make sure that we're talking about the same fliptime
firstFlipDiff = fwdFlipTimesR(firstFwdIndR)-fwdFlipTimesL(firstFwdIndL) ;
if firstFlipDiff > .001
    firstFwdIndL = firstFwdIndL + 1;
elseif firstFlipDiff < -.001
    firstFwdIndR = firstFwdIndR + 1;
end

%how many strokes to look at after perturbation
maxWingstrokes = min([(length(fwdFlipTimesR(firstFwdIndR:end))-1),...
    (length(fwdFlipTimesL(firstFwdIndL:end))-1)]) ;

flipTimesR = fwdFlipTimesR(firstFwdIndR:(firstFwdIndR+maxWingstrokes));
flipTimesL = fwdFlipTimesL(firstFwdIndL:(firstFwdIndL+maxWingstrokes));
flipTimes = mean([flipTimesR ;flipTimesL],1) ;

phi_frontR = fnval(sp_phiR, flipTimesR) ; 
phi_frontL = fnval(sp_phiL, flipTimesL) ; 

Phi_Front = mean([phi_frontR; phi_frontL],1) ;

%% Find average values of phi_front to compute the difference

prePertIndR = find((fwdFlipTimesR > correctionTime(1)/1000) & (fwdFlipTimesR < 0)) ;
prePertIndL = find((fwdFlipTimesL > correctionTime(1)/1000) & (fwdFlipTimesL < 0)) ;

prePertTimesR = fwdFlipTimesR(prePertIndR) ;
prePertTimesL = fwdFlipTimesL(prePertIndL) ;

phi_frontR_avg = mean(fnval(sp_phiR, prePertTimesR)) ; 
phi_frontL_avg = mean(fnval(sp_phiL, prePertTimesL)) ; 

Phi_FrontAvg = mean([phi_frontR_avg, phi_frontL_avg]) ; 

Phi_FrontDiff = Phi_Front - Phi_FrontAvg ; 

%% Find error in control predictions

parfor q = 1:systemSize
    [i,j,k,m] = ind2sub(dim,q) ;
    controlPred = K(i,j,k,m) + Kp(i,j,k,m)*fnval(sp_velocity,flipTimes-DeltaT(i,j,k,m))...
        + Ki(i,j,k,m)*(fnval(sp_pitch,flipTimes-DeltaT(i,j,k,m)) - bodyPitchAvg) ;
    
    Err = Phi_FrontDiff - controlPred;
    MSErr = mean(Err.^2) ;
    
    ErrorList(q) = MSErr ;
end

ErrorCube = reshape(ErrorList,dim) ; %I can't not use this name, right?

[MinErr, MinErrInd] = min(ErrorList) ; 
[Idx1, Idx2, Idx3, Idx4] = ind2sub(dim, MinErrInd) ;
%MinErr
minKi = Ki(Idx1, Idx2, Idx3, Idx4); 
minKp = Kp(Idx1, Idx2, Idx3, Idx4); 
minDt = DeltaT(Idx1, Idx2, Idx3, Idx4);
%minDt2 = DeltaT2(Idx1, Idx2, Idx3, Idx4, Idx5);
minK = K(Idx1, Idx2, Idx3, Idx4); 
results_full = [minKi, minKp, minDt, minK, MinErr] ;

%% Check other control models
%Not adding a constant term
%{
NoKInd = find(K == 0) ;
[MinErrNoK, MinErrIndNoK] = min(ErrorList(NoKInd)) ;
MinErrIndNoK = MinErrIndNoK + NoKInd(1) - 1 ;
[Idx1, Idx2, Idx3, Idx4, Idx5] = ind2sub(dim, MinErrIndNoK) ;
minKi_NoK = Ki(Idx1, Idx2, Idx3, Idx4, Idx5); 
minKp_NoK = Kp(Idx1, Idx2, Idx3, Idx4, Idx5); 
minDt_NoK = DeltaT(Idx1, Idx2, Idx3, Idx4, Idx5);
minDt2_NoK = DeltaT2(Idx1, Idx2, Idx3, Idx4, Idx5);
minK_NoK = K(Idx1, Idx2, Idx3, Idx4, Idx5); 
results_NoK = [minKi_NoK, minKp_NoK, minDt_NoK, minDt2_NoK, minK_NoK, MinErrNoK] ;

%No integral term
NoKiInd = find(Ki == 0) ;
[MinErrNoKi, MinErrIndNoKi] = min(ErrorList(NoKiInd)) ;
MinErrIndNoKi = ki_N*(MinErrIndNoKi(1) - 1) + 1 ;
[Idx1, Idx2, Idx3, Idx4, Idx5] = ind2sub(dim, MinErrIndNoKi) ;
minKi_NoKi = Ki(Idx1, Idx2, Idx3, Idx4, Idx5); 
minKp_NoKi = Kp(Idx1, Idx2, Idx3, Idx4, Idx5); 
minDt_NoKi = DeltaT(Idx1, Idx2, Idx3, Idx4, Idx5);
minDt2_NoKi = DeltaT2(Idx1, Idx2, Idx3, Idx4, Idx5);
minK_NoKi = K(Idx1, Idx2, Idx3, Idx4, Idx5); 
results_NoKi = [minKi_NoKi, minKp_NoKi, minDt_NoKi, minDt2_NoKi, minK_NoKi, MinErrNoKi] ;

%Fix the time delays to be the same for the P and I terms
SameTInd = find(DeltaT == DeltaT2) ; 
temp = nan(size(ErrorList)) ; 
temp(SameTInd) = ErrorList(SameTInd) ; 
[MinErrSameT, MinErrIndSameT] = min(temp) ;
[Idx1, Idx2, Idx3, Idx4, Idx5] = ind2sub(dim, MinErrIndSameT) ;
minKi_SameT = Ki(Idx1, Idx2, Idx3, Idx4, Idx5); 
minKp_SameT = Kp(Idx1, Idx2, Idx3, Idx4, Idx5); 
minDt_SameT = DeltaT(Idx1, Idx2, Idx3, Idx4, Idx5);
minDt2_SameT = DeltaT2(Idx1, Idx2, Idx3, Idx4, Idx5);
minK_SameT = K(Idx1, Idx2, Idx3, Idx4, Idx5); 
results_SameT = [minKi_SameT, minKp_SameT, minDt_SameT, minDt2_SameT, minK_SameT, MinErrSameT] ;

%Fix the time delays to be the same for the P and I terms and no constant
SameTNoKInd = intersect(SameTInd,NoKInd) ; 
temp = nan(size(ErrorList)) ; 
temp(SameTNoKInd) = ErrorList(SameTNoKInd) ; 
[MinErrSameTNoK, MinErrIndSameTNoK] = min(temp) ;
[Idx1, Idx2, Idx3, Idx4, Idx5] = ind2sub(dim, MinErrIndSameTNoK) ;
minKi_SameT_NoK = Ki(Idx1, Idx2, Idx3, Idx4, Idx5); 
minKp_SameT_NoK = Kp(Idx1, Idx2, Idx3, Idx4, Idx5); 
minDt_SameT_NoK = DeltaT(Idx1, Idx2, Idx3, Idx4, Idx5);
minDt2_SameT_NoK = DeltaT2(Idx1, Idx2, Idx3, Idx4, Idx5);
minK_SameT_NoK = K(Idx1, Idx2, Idx3, Idx4, Idx5); 
results_SameT_NoK = [minKi_SameT_NoK, minKp_SameT_NoK, minDt_SameT_NoK, minDt2_SameT_NoK, minK_SameT_NoK, MinErrSameTNoK] ;

%No proportional term
NoKpInd = find(Kp == 0) ; 
temp = nan(size(ErrorList)) ; 
temp(NoKpInd) = ErrorList(NoKpInd) ; 
[MinErrNoKp, MinErrIndNoKpInd] = min(temp) ;
[Idx1, Idx2, Idx3, Idx4, Idx5] = ind2sub(dim, MinErrIndNoKpInd) ;
minKi_NoKp = Ki(Idx1, Idx2, Idx3, Idx4, Idx5); 
minKp_NoKp = Kp(Idx1, Idx2, Idx3, Idx4, Idx5); 
minDt_NoKp = DeltaT(Idx1, Idx2, Idx3, Idx4, Idx5);
minDt2_NoKp = DeltaT2(Idx1, Idx2, Idx3, Idx4, Idx5);
minK_NoKp = K(Idx1, Idx2, Idx3, Idx4, Idx5); 
results_NoKp = [minKi_NoKp, minKp_NoKp, minDt_NoKp, minDt2_NoKp, minK_NoKp, MinErrNoKp] ;

%No proportional term, no constant
NoKpNoKInd = intersect(NoKpInd,NoKInd) ; 
temp = nan(size(ErrorList)) ; 
temp(NoKpNoKInd) = ErrorList(NoKpNoKInd) ; 
[MinErrNoKpNoK, MinErrIndNoKpNoKInd] = min(temp) ;
[Idx1, Idx2, Idx3, Idx4, Idx5] = ind2sub(dim, MinErrIndNoKpNoKInd) ;
minKi_NoKp_NoK = Ki(Idx1, Idx2, Idx3, Idx4, Idx5); 
minKp_NoKp_NoK = Kp(Idx1, Idx2, Idx3, Idx4, Idx5); 
minDt_NoKp_NoK = DeltaT(Idx1, Idx2, Idx3, Idx4, Idx5);
minDt2_NoKp_NoK = DeltaT2(Idx1, Idx2, Idx3, Idx4, Idx5);
minK_NoKp_NoK = K(Idx1, Idx2, Idx3, Idx4, Idx5); 
results_NoKp_NoK = [minKi_NoKp_NoK, minKp_NoKp_NoK, minDt_NoKp_NoK, minDt2_NoKp_NoK, minK_NoKp_NoK, MinErrNoKpNoK] ;
%}
%% Plot slices of controller space

if plotFlag2
    %{
    for s = 1:10:deltat_N
       figure;
       C = contour(Ki(:,:,s),Kp(:,:,s),ErrorCube(:,:,s)) ;
       xlabel('K_i')
       ylabel('K_p')
       title(['\Delta t = ' num2str(1000*deltat_lin(s),4) ' [ms]'])
       clabel(C) ; 
    end
    %}
    h1 = figure ;
    C_fixT = contour(Ki(:,:,Idx3,Idx4),Kp(:,:,Idx3,Idx4),ErrorCube(:,:,Idx3,Idx4)) ;
    hold on
    [C_Err_fixT, h_err_fixT] = contour(Ki(:,:,Idx3,Idx4),Kp(:,:,Idx3,Idx4),ErrorCube(:,:,Idx3,Idx4),[MinErr+4 MinErr+4]) ; 
    xlabel('K_i')
    ylabel('K_p')
    title(['\Delta t = ' num2str(1000*minDt,4) ' [ms], K = ' num2str(minK)])
    clabel(C_fixT) ; 
    
    Ki_Dist1 = abs(C_Err_fixT(1,2:end) - minKi) ; 
    Kp_Dist1 = abs(C_Err_fixT(2,2:end) - minKp) ; 
    
    h2 = figure ;
    C_fixKi = contour(squeeze(DeltaT(Idx1,:,:,Idx4)),squeeze(Kp(Idx1,:,:,Idx4)),squeeze(ErrorCube(Idx1,:,:,Idx4))) ;
    hold on
    [C_Err_fixKi, h_err_fixKi] = contour(squeeze(DeltaT(Idx1,:,:,Idx4)),squeeze(Kp(Idx1,:,:,Idx4)),squeeze(ErrorCube(Idx1,:,:,Idx4)),[MinErr+4 MinErr+4]) ; 
    xlabel('\Delta T')
    ylabel('K_p')
    title(['K_i = ' num2str(minKi,4) ' , K = ' num2str(minK)])
    clabel(C_fixKi) ; 
    
    Dt_Dist1 = abs(C_Err_fixKi(1,2:end) - minDt) ; 
    Kp_Dist2 = abs(C_Err_fixKi(2,2:end) - minKp) ; 
    
    h3 = figure ;
    C_fixKp = contour(squeeze(Ki(:,Idx2,:,Idx4)),squeeze(DeltaT(:,Idx2,:,Idx4)),squeeze(ErrorCube(:,Idx2,:,Idx4))) ;
    hold on
    [C_Err_fixKp, h_err_fixKp] = contour(squeeze(Ki(:,Idx2,:,Idx4)),squeeze(DeltaT(:,Idx2,:,Idx4)),squeeze(ErrorCube(:,Idx2,:,Idx4)),[MinErr+4 MinErr+4]) ; 
    xlabel('K_i')
    ylabel('\Delta T')
    title(['K_p = ' num2str(1000*minKp,4) ' [ms] , K = ' num2str(minK)])
    clabel(C_fixKp) ; 
    
    Ki_Dist2 = abs(C_Err_fixKp(1,2:end) - minKi) ; 
    Dt_Dist2 = abs(C_Err_fixKp(2,2:end) - minDt) ;
    
    Ki_var = max([Ki_Dist1, Ki_Dist2]) ;
    Kp_var = max([Kp_Dist1, Kp_Dist2]) ;
    Dt_var = max([Dt_Dist1, Dt_Dist2]) ;
    
    CI = [Ki_var, Kp_var, Dt_var] ;
    
    %close(h1) ; close(h2) ; close(h3) ;
end


%% Plot controller fit
if plotFlag3
    %{
    hctrlfit = figure ;
    set(hctrlfit,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' CtrlFit'],'numbertitle','off')
    hold on
    controlPred = minK + minKp*fnval(sp_velocity,t_AC-minDt)...
        + minKi*(fnval(sp_pitch,t_AC-minDt2)- bodyPitchAvg) ; 
    plot(t_AC*1000, controlPred + Phi_FrontAvg, 'r-','LineWidth',2) ; 
    plot(flipTimes*1000, Phi_Front, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
    set(gca,'xlim',[0 max(flipTimes)*1000+10])
    legend({'Controller','Data'},'location','northeast')
    xlabel('Time [ms]')
    ylabel('Front Stroke Angle [deg]')
    title(['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum)])
    %}
    controlPred = minK + minKp*fnval(sp_velocity,t_AC-minDt)...
        + minKi*(fnval(sp_pitch,t_AC-minDt)- bodyPitchAvg) ;
    %controlPredErr = maxK + (maxKp)*fnval(sp_velocity,t_AC-(maxDt))...
    %    + (maxKi)*(fnval(sp_pitch,(maxDt))- bodyPitchAvg) ;
    %{
    controlPred_NoK = minKp_NoK*fnval(sp_velocity,t_AC-minDt_NoK)...
        + minKi_NoK*(fnval(sp_pitch,t_AC-minDt2_NoK)- bodyPitchAvg) ; 
    
    controlPred_NoKi = minK_NoKi + minKp_NoKi*fnval(sp_velocity,t_AC-minDt_NoKi) ;
    
    controlPred_NoKp = minK_NoKp + minKi_NoKp*(fnval(sp_pitch,t_AC-minDt2_NoKp)- bodyPitchAvg) ; 
    
    controlPred_NoKp_NoK = minKi_NoKp_NoK*(fnval(sp_pitch,t_AC-minDt2_NoKp_NoK)- bodyPitchAvg) ; 
    
    controlPred_SameT = minK_SameT + minKp_SameT*fnval(sp_velocity,t_AC-minDt_SameT)...
        + minKi_SameT*(fnval(sp_pitch,t_AC-minDt2_SameT)- bodyPitchAvg) ; 
    
    controlPred_SameT_NoK = minKp_SameT_NoK*fnval(sp_velocity,t_AC-minDt_SameT_NoK)...
        + minKi_SameT_NoK*(fnval(sp_pitch,t_AC-minDt2_SameT_NoK)- bodyPitchAvg) ; 
    %}
    %xlim = [1000*flipTimes(1)-2 1000*flipTimes(end)+2] ; 
    ylim = [-20 20] ;
    tsfvec = [0 7 7 0 0] ;
    
    hcontrib = figure ;
    hold on
    %set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'Position', [500 500 420 140]);
    set(gcf, 'Color', 'w')
    set(gcf, 'Renderer','opengl')
    set(hcontrib,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' Contributions'],'numbertitle','off')
    
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',.7) ;
    set(hf,'HandleVisibility','off')
   
    hP =  plot(t_AC*1000,minKp*fnval(sp_velocity,t_AC - minDt),'LineStyle','-','Color',[1 1 1]*.5,'LineWidth',1);
    hI = plot(t_AC*1000,minKi*(fnval(sp_pitch,t_AC - minDt)- bodyPitchAvg),'LineStyle','--','Color',[101 67 33]/256,'LineWidth',1);
    shadedErrorBar(t_AC*1000, controlPred, 2*ones(size(controlPred)),{'-','LineWidth',2.5,'Color',[96 96 255]/256},1) ;
    %hPI = plot(t_AC*1000, controlPred, '-','LineWidth',2,'Color',[96 96 255]/256) ;    
    hdata = plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', [1 .5 0 ], 'MarkerSize',6);
    %plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
     
    %set(gca,'xlim',[-10 43.5])
    %set(gca,'ylim',ylim)
    set(gca,'fontsize',12)
    %savefig('controller')
    %xlabel('Time [ms]','fontsize',10)
    %ylabel({'Front Stroke Angle Deviation';'  [deg]'},'fontsize',10)
    %legend([hdata,hPI, hP, hI],{'\Delta\phi_{data}','\Delta\phi_{ctrl}','K_p\omega(t-\DeltaT)',...
    %    'K_i\Delta\theta(t-\DeltaT)'},'location','northeast')
    %title('PI Control with Constant')
    %{
    subplot(3,1,1)
        plot(t_AC*1000,fnval(sp_pitch,t_AC)- bodyPitchAvg,'k.')
        hold on
        plot(t_AC*1000,fnval(sp_pitch,t_AC - minDt2)- bodyPitchAvg,'rx') 
        xlim = get(gca,'xlim') ;
        xlabel('Time [ms]')
        ylabel('Pitch Angle Deviation [deg]')
        legend({'\Delta\theta(t)','\Delta\theta(t-\DeltaT_2)'},'location','northeast')
        title(['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum)])
    
    subplot(3,1,2)
    %subplot(2,1,1)
        plot(t_AC*1000,fnval(sp_velocity,t_AC),'k.')
        hold on
        plot(t_AC*1000,fnval(sp_velocity,t_AC - minDt),'rx')
        set(gca,'xlim',xlim)
        xlabel('Time [ms]')
        ylabel('Pitch Velocity [deg/s]')
        legend({'\omega(t)','\omega(t-\DeltaT_1)'},'location','northeast')
        %title('Pitch and delayed pitch')
    subplot(2,2,1)
     %subplot(2,1,2)
        plot(t_AC*1000, controlPred_NoK, 'r-','LineWidth',2.5)
        hold on
        plot(t_AC*1000,minKp_NoK*fnval(sp_velocity,t_AC-minDt_NoK),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        plot(t_AC*1000,minKi_NoK*(fnval(sp_pitch,t_AC-minDt2_NoK)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}','K_p\omega(t-\DeltaT_1)','K_i\Delta\theta(t-\DeltaT_2)'...
            '\Delta\phi_{data}'},'location','northeast')
        title('No Constant')
    subplot(2,2,2)
     %subplot(2,1,2)
        plot(t_AC*1000, controlPred_NoKi, 'r-','LineWidth',2.5)
        hold on
        plot(t_AC*1000,minKp_NoKi*fnval(sp_velocity,t_AC-minDt_NoKi),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        %plot(t_AC*1000,minKi_NoK*(fnval(sp_pitch,t_AC-minDt2_NoK)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}','K_p\omega(t-\DeltaT_1)',...%'K_i\Delta\theta(t-\DeltaT_2)'...
            '\Delta\phi_{data}'},'location','northeast')
        title('P Control')
    subplot(2,2,[3,4])
     %subplot(2,1,2)
        plot(t_AC*1000, controlPred, 'r-','LineWidth',2.5)
        hold on
        plot(t_AC*1000,minKp*fnval(sp_velocity,t_AC - minDt),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        plot(t_AC*1000,minKi*(fnval(sp_pitch,t_AC - minDt2)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}','K_p\omega(t-\DeltaT_1)','K_i\Delta\theta(t-\DeltaT_2)'...
            '\Delta\phi_{data}'},'location','northeast')
        title('PI Control with Constant')
    
        
    hcontrib2 = figure ;
    set(hcontrib2,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' Contributions 2'],'numbertitle','off')
    
    subplot(2,2,1)
        plot(t_AC*1000, controlPred_SameT_NoK, 'r-','LineWidth',2.5)
        hold on
        plot(t_AC*1000,minKp_SameT_NoK*fnval(sp_velocity,t_AC-minDt_SameT_NoK),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        plot(t_AC*1000,minKi_SameT_NoK*(fnval(sp_pitch,t_AC-minDt2_SameT_NoK)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}','K_p\omega(t-\DeltaT)','K_i\Delta\theta(t-\DeltaT)'...
            '\Delta\phi_{data}'},'location','northeast')
        title('Same Time Delay, No Constant')
    subplot(2,2,2)
        plot(t_AC*1000, controlPred_NoKp, 'r-','LineWidth',2.5)
        hold on
        %plot(t_AC*1000,minKp_NoKi*fnval(sp_velocity,t_AC-minDt_NoKi),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        plot(t_AC*1000,minKi_NoKp*(fnval(sp_pitch,t_AC-minDt2_NoKp)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}',...%'K_p\omega(t-\DeltaT_1)',...%
            'K_i\Delta\theta(t-\DeltaT_2)','\Delta\phi_{data}'},'location','northeast')
        title('I Control with Constant')
    subplot(2,2,3)
        plot(t_AC*1000, controlPred_SameT, 'r-','LineWidth',2.5)
        hold on
        plot(t_AC*1000,minKp_SameT*fnval(sp_velocity,t_AC - minDt_SameT),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        plot(t_AC*1000,minKi_SameT*(fnval(sp_pitch,t_AC - minDt2_SameT)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}','K_p\omega(t-\DeltaT)','K_i\Delta\theta(t-\DeltaT)'...
            '\Delta\phi_{data}'},'location','northeast')
        title('Same Time Delay with Constant')
    subplot(2,2,4)
        plot(t_AC*1000, controlPred_NoKp_NoK, 'r-','LineWidth',2.5)
        hold on
        %plot(t_AC*1000,minKp*fnval(sp_velocity,t_AC - minDt),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        plot(t_AC*1000,minKi_NoKp_NoK*(fnval(sp_pitch,t_AC - minDt2_NoKp_NoK)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}',...%'K_p\omega(t-\DeltaT_1)',
            'K_i\Delta\theta(t-\DeltaT_2)'...
            '\Delta\phi_{data}'},'location','northeast')
        title('I Control without Constant')
    %}
end   



