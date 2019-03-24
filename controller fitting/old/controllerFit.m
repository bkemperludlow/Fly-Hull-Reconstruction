function [minKi, minKp, minDt, MinErr, PitchEstErr] = ...
    controllerFit(ExprNum,MovNum,PitchType,varyPitchSplineFlag,...
        plotFlag,plotFlag2,plotFlag3)
%Create 3D space of control parameters, find minimum

defineConstantsScript
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
    datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
        num2str(ExprNum), '_mov_',zstr,num2str(MovNum)) ;
elseif strcmp(PitchType,'down')
    datapath = strcat('F:\luca\Analysis\pitch down\Expr_',...
        num2str(ExprNum), '_mov_',zstr,num2str(MovNum)) ;
else
    disp('check PitchType')
    return ;
end
cd(datapath)

datafilename = strcat(datapath,'\Expr',num2str(ExprNum), ...
    'mov',zstr,num2str(MovNum), '_Data_manually_corrected.mat') ;

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
    pulseLength = 8 ; %ms
end

%Get flip times
if isfield(data,'fwdFlipTimesR') && isfield(data,'backFlipTimesR')
    fwdFlipTimesR = data.fwdFlipTimesR ;
    %backFlipTimesR = data.backFlipTimesR ;
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
    maxRMS = .5 ;
    
    [sp_pitch, pitch_smooth, RMS, PitchEstErr] = varySpline(currtvec,currBodyPitch,...
        PitchErrLow,PitchErrHigh,Ntest,maxRMS,plotFlag); 
else
    PitchEstErr = .4 ;
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

ki_min = -20 ; %unitless
ki_max = 20 ;
ki_N = 50 ;

kp_min = 0 ; %seconds
kp_max = .1 ; 
kp_N = 100 ;

deltat_min = .004 ; %seconds
deltat_max = .03 ;
deltat_N = 100 ; 

ki_lin = linspace(ki_min, ki_max, ki_N) ;
kp_lin = linspace(kp_min, kp_max, kp_N) ;
deltat_lin = linspace(deltat_min, deltat_max, deltat_N) ;

systemSize = ki_N*kp_N*deltat_N ;

[Ki,Kp,DeltaT] = meshgrid(ki_lin,kp_lin,deltat_lin) ; 
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
    [i,j,k] = ind2sub(dim,q) ;
    controlPred = Kp(i,j,k)*fnval(sp_velocity,flipTimes-DeltaT(i,j,k))...
        + Ki(i,j,k)*(fnval(sp_pitch,flipTimes-DeltaT(i,j,k)) - bodyPitchAvg) ;
    
    Err = Phi_FrontDiff - controlPred;
    RMSErr = sqrt(mean(Err.^2)) ;
    
    ErrorList(q) = RMSErr ;
end

ErrorCube = reshape(ErrorList,dim) ; %I can't not use this name, right?

[MinErr, MinErrInd] = min(ErrorList) ; 
[Idx1, Idx2, Idx3] = ind2sub(dim, MinErrInd) ;
%MinErr
minKi = Ki(Idx1, Idx2, Idx3); 
minKp = Kp(Idx1, Idx2, Idx3); 
minDt = DeltaT(Idx1, Idx2, Idx3);

%% Plot slices of controller space

if plotFlag2
    for s = 1:10:deltat_N
       figure;
       C = contour(Ki(:,:,s),Kp(:,:,s),ErrorCube(:,:,s)) ;
       xlabel('K_i')
       ylabel('K_p')
       title(['\Delta t = ' num2str(1000*deltat_lin(s),4) ' [ms]'])
       clabel(C) ; 
    end
    
    figure ;
    C = contour(Ki(:,:,Idx3),Kp(:,:,Idx3),ErrorCube(:,:,Idx3)) ;
    xlabel('K_i')
    ylabel('K_p')
    title(['Optimal \Delta t = ' num2str(1000*deltat_lin(Idx3),4) ' [ms]'])
    clabel(C) ; 
end

if plotFlag3
    hctrlfit = figure ;
    set(hctrlfit,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' CtrlFit'],'numbertitle','off')
    hold on
    controlPred = minKp*fnval(sp_velocity,t_AC-minDt)...
        + minKi*(fnval(sp_pitch,t_AC-minDt)- bodyPitchAvg) ; 
    plot(t_AC*1000, controlPred + Phi_FrontAvg, 'r-','LineWidth',2) ; 
    plot(flipTimes*1000, Phi_Front, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
    set(gca,'xlim',[0 max(flipTimes)*1000+10])
    legend({'Controller','Data'},'location','northeast')
    xlabel('Time [ms]')
    ylabel('Front Stroke Angle [deg]')
    title(['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum)])
    
    hcontrib = figure ;
    set(hcontrib,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' Contributions'],'numbertitle','off')
    
    subplot(3,1,1)
        plot(t_AC*1000,fnval(sp_pitch,t_AC)- bodyPitchAvg,'k.')
        hold on
        plot(t_AC*1000,fnval(sp_pitch,t_AC - minDt)- bodyPitchAvg,'rx') 
        xlim = get(gca,'xlim') ;
        xlabel('Time [ms]')
        ylabel('Pitch Angle Deviation [deg]')
        legend({'\Delta\theta(t)','\Delta\theta(t- \Delta T'},'location','northeast')
        title(['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum)])
    
    subplot(3,1,2)
    %subplot(2,1,1)
        plot(t_AC*1000,fnval(sp_velocity,t_AC),'k.')
        hold on
        plot(t_AC*1000,fnval(sp_velocity,t_AC - minDt),'rx')
        set(gca,'xlim',xlim)
        xlabel('Time [ms]')
        ylabel('Pitch Velocity [deg/s]')
        legend({'\omega(t)','\omega(t- \Delta T'},'location','northeast')
        %title('Pitch and delayed pitch')
     subplot(3,1,3)
     %subplot(2,1,2)
        plot(t_AC*1000, controlPred, 'r-','LineWidth',2.5)
        hold on
        plot(t_AC*1000,minKp*fnval(sp_velocity,t_AC - minDt),'LineStyle','--','Color',[0 .9 .9],'LineWidth',1.5)
        plot(t_AC*1000,minKi*(fnval(sp_pitch,t_AC - minDt)- bodyPitchAvg) ,'m--','LineWidth',1.5)
        plot(flipTimes*1000, Phi_Front - Phi_FrontAvg, 'ko', 'MarkerFaceColor', 'g', 'MarkerSize',10)
        set(gca,'xlim',xlim)
        set(gca,'fontsize',12)
        xlabel('Time [ms]','fontsize',14)
        ylabel('Front Stroke Angle Deviation [deg]','fontsize',14)
        legend({'\Delta\phi_{ctrl}','K_p\omega(t-\DeltaT)',...%'K_i\Delta\theta(t-\DeltaT)'...
            '\Delta\phi_{data}'},'location','northeast')
end   



