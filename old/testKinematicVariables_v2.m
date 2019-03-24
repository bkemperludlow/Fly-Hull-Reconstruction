%% Preliminaries
cd('F:\luca\Analysis\pitch down\Expr_7_mov_040')
load('Expr7mov040_Data_manually_corrected.mat')
defineConstantsScript
pitchType = 'down' ;

phiCut = false ;
etaCut = false;
thetaCut = true ;

savePath = 'C:\Users\Fruit Flies\Desktop\torque questions\Cut and Paste Kinematics\pitch down (7.40)' ;
saveStr = ['PhiCut = ' num2str(phiCut) ', EtaCut = ' num2str(etaCut) ', ThetaCut = ' num2str(thetaCut) ] ;
saveFlag = false ;

EstPhiErr = 1 ;
EstEtaErr = 2.5 ; %2.5
EstThetaErr = 2 ; %2

checkPlotFlag1 = true ;
checkPlotFlag2 = true ;
%------------------------------------------------------------------------

if isfield(data,'t')
    t = data.t ;
else
    t1_temp = data.params.startTrackingTime ;
    t2_temp = data.params.endTrackingTime ;
    fps = data.params.fps ;
    t = (t1_temp:t2_temp)/fps ;
end
tms = t*1000 ;
Nframes = length(t) ;

if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

if isfield(data, 'manualCorrRangeMS')
    xlim = data.manualCorrRangeMS ;
elseif isfield(data, 'correctionTime')
    xlim = data.correctionTime ;
else
    xlim = [tms(1) tms(end)] ;
end

%% Read relevant data
backFlipTimesR = data.backFlipTimesR ;
fwdFlipTimesR = data.fwdFlipTimesR ;
%backFlipTimesL = data.backFlipTimesL ;
%fwdFlipTimesL = data.fwdFlipTimesL ;

phiR = -data.anglesBodyFrame(:,PHIR) ;
etaR = (180/pi)*unwrap((pi/180)*data.anglesBodyFrame(:,ETAR)) ;
thetaR = data.anglesBodyFrame(:,THETAR) ;
phiR(ignoreFrames) = NaN ; etaR(ignoreFrames) = NaN ; thetaR(ignoreFrames) = NaN ;

phiL = data.anglesBodyFrame(:,PHIL) ;
etaL = (180/pi)*unwrap((pi/180)*data.anglesBodyFrame(:,ETAL)) ;
thetaL = data.anglesBodyFrame(:,THETAL) ;
phiL(ignoreFrames) = NaN ; etaL(ignoreFrames) = NaN ; thetaL(ignoreFrames) = NaN ;

%bodyCM = data.bodyCM ;
thetaB = data.anglesLabFrame(:,THETAB) ;
%[maxPitch, PitchInd] = max(thetaB) ;
%PitchInd = PitchInd ; %- 40 ;
[minPitch, PitchInd] = min(thetaB) ;

FlipInd1 = find(backFlipTimesR < t(PitchInd),1,'last') ;
FlipInd2 = find(backFlipTimesR > t(PitchInd),1,'first') ;
delT = backFlipTimesR(FlipInd2) - backFlipTimesR(FlipInd1) ; 
Ind1 = find(t == backFlipTimesR(FlipInd1)) ;
Ind2 = find(t == backFlipTimesR(FlipInd2)) ;

Nframes = 50 ;
repNum = 5 ;

currt = linspace(0, delT*repNum, Nframes*repNum) ;
%currt = t(Ind1:Ind2) ;

%% Find times interval to cut and paste (cut [cutpoint:zeropoint] and paste afterwards)
prePertFlipInds = find(backFlipTimesR < 0,2,'last') ;
cutInd1 = find(t == backFlipTimesR(prePertFlipInds(1))) ;
cutInd2 = find(t == backFlipTimesR(prePertFlipInds(2))) ;

if phiCut
    currPhiR = repmat(phiR(cutInd1:cutInd2),repNum,1) ;
    currPhiL = repmat(phiL(cutInd1:cutInd2),repNum,1) ;
else
    [PhiInd1, PhiInd2] = chooseCutWindow(phiR,phiL,t,[Ind1 Ind2],...
        backFlipTimesR,fwdFlipTimesR) ;
    currPhiR = repmat(phiR(PhiInd1:PhiInd2),repNum,1) ;
    currPhiL = repmat(phiL(PhiInd1:PhiInd2),repNum,1) ;
end

if etaCut
    currEtaR = repmat(etaR(cutInd1:cutInd2),repNum,1) ;
    currEtaL = repmat(etaL(cutInd1:cutInd2),repNum,1) ;
else
    [EtaInd1, EtaInd2] = chooseCutWindow(etaR,etaL,t,[Ind1 Ind2],...
        backFlipTimesR,fwdFlipTimesR) ;
    
    currEtaR = repmat(etaR(EtaInd1:EtaInd2),repNum,1) ;
    currEtaL = repmat(etaL(EtaInd1:EtaInd2),repNum,1) ;
end

if thetaCut
    currThetaR = repmat(thetaR(cutInd1:cutInd2),repNum,1) ;
    currThetaL = repmat(thetaL(cutInd1:cutInd2),repNum,1) ;
else
    [ThetaInd1, ThetaInd2] = chooseCutWindow(thetaR,thetaL,t,[Ind1 Ind2],...
        backFlipTimesR,fwdFlipTimesR) ;
    
    currThetaR = repmat(thetaR(ThetaInd1:ThetaInd2),repNum,1) ;
    currThetaL = repmat(thetaL(ThetaInd1:ThetaInd2),repNum,1) ;
end

%% Smooth wing angles
currtvecPhi = linspace(0, delT*repNum, length(currPhiR)) ;
[sp_phiR, ~, ~] = mySplineSmooth(currtvecPhi,currPhiR,EstPhiErr) ;
[sp_phiL, ~, ~] = mySplineSmooth(currtvecPhi,currPhiL,EstPhiErr) ;
phiR_smooth = fnval(sp_phiR,currt) ;
phiL_smooth = fnval(sp_phiL,currt) ;

currtvecEta = linspace(0, delT*repNum, length(currEtaR)) ;
[sp_etaR, ~, ~] = mySplineSmooth(currtvecEta,currEtaR,EstEtaErr) ;
[sp_etaL, ~, ~] = mySplineSmooth(currtvecEta,currEtaL,EstEtaErr) ;
etaR_smooth = fnval(sp_etaR,currt) ;
etaL_smooth = fnval(sp_etaL,currt) ;

currtvecTheta = linspace(0, delT*repNum, length(currThetaR)) ;
[sp_thetaR, ~, ~] = mySplineSmooth(currtvecTheta,currThetaR,EstThetaErr) ;
[sp_thetaL, ~, ~] = mySplineSmooth(currtvecTheta,currThetaL,EstThetaErr) ;
thetaR_smooth = fnval(sp_thetaR,currt) ;
thetaL_smooth = fnval(sp_thetaL,currt) ;

if checkPlotFlag1
    hkinematics = figure ;
    subplot(3,1,1)
        hold on
        plot(currt*1000,phiR_smooth,'r-')
        plot(currt*1000,phiL_smooth,'b-')
        ylabel('\phi [deg]')
    subplot(3,1,2)
        hold on
        plot(currt*1000,etaR_smooth,'r-')
        plot(currt*1000,etaL_smooth,'b-')
        ylabel('\eta [deg]')
    subplot(3,1,3)
        hold on
        plot(currt*1000,thetaR_smooth,'r-')
        plot(currt*1000,thetaL_smooth,'b-')
        ylabel('\theta [deg]')
        xlabel('Time [ms]')
end
%{
if phiCut
    newFlipTimesR =  mean(fnzeros(fnder(sp_phiR,1), [t(1) t(end)]),1) ;
    %figure ; plot(t,phiR_smooth) ; hold on ; 
    %plot(phiRdotZeros,180*ones(size(newFlipTimes)), 'ro')
    if fnval(sp_phiR,newFlipTimesR(1)) > 100
        data.backFlipTimesR = newFlipTimesR(1:2:end) ;
        data.fwdFlipTimesR = newFlipTimesR(2:2:end) ;
    else
        data.backFlipTimesR = newFlipTimesR(2:2:end) ;
        data.fwdFlipTimesR = newFlipTimesR(1:2:end) ;
    end
    
    newFlipTimesL =  fnzeros(fnder(sp_phiL,1), [t(1) t(end)]) ;
    %figure ; plot(t,phiR_smooth) ; hold on ; 
    %plot(phiRdotZeros,180*ones(size(newFlipTimes)), 'ro')
    if fnval(sp_phiR,newFlipTimesL(1)) > 100
        data.backFlipTimesL = newFlipTimesL(1:2:end) ;
        data.fwdFlipTimesL = newFlipTimesL(2:2:end) ;
    else
        data.backFlipTimesL = newFlipTimesL(2:2:end) ;
        data.fwdFlipTimesL = newFlipTimesL(1:2:end) ;
    end
        
end
%}
%% Define wing vectors based on Euler angles

%these are the vector hats for the wings assuming that the stroke plane is
%coincident with the xy plane and the body axis is pointing along the y
%direction
spanHatR_xy = [sind(phiR_smooth).*cosd(thetaR_smooth); cosd(phiR_smooth).*cosd(thetaR_smooth); sind(thetaR_smooth)] ;
spanHatL_xy = [-sind(phiL_smooth).*cosd(thetaL_smooth); cosd(phiL_smooth).*cosd(thetaL_smooth); sind(thetaL_smooth)] ;

chordHatR_xy = [-cosd(etaR_smooth).*cosd(phiR_smooth); cosd(etaR_smooth).*sind(phiR_smooth) ; sind(etaR_smooth) ] ;
chordHatL_xy = [cosd(etaL_smooth).*cosd(phiL_smooth); cosd(etaL_smooth).*sind(phiL_smooth) ; sind(etaL_smooth) ] ;

%{
badIndR = find(spanHatR_xy(1,:) < 0) ; 
badIndL = find(spanHatL_xy(1,:) > 0) ;
badInd = unique([badIndR, badIndL]) ;

NanMat = repmat([NaN; NaN; NaN], 1, length(badInd)) ;
spanHatR_xy(:,badInd) = NanMat ;
spanHatL_xy(:,badInd) = NanMat ;
chordHatR_xy(:,badInd) = NanMat ;
chordHatL_xy(:,badInd) = NanMat ;
%}
if checkPlotFlag2
    hspans = figure; 
    plot3(spanHatR_xy(1,:),spanHatR_xy(2,:),spanHatR_xy(3,:),'ko','MarkerFaceColor','r')
    hold on
    plot3(spanHatL_xy(1,:),spanHatL_xy(2,:),spanHatL_xy(3,:),'ko','MarkerFaceColor','b')
    axis equal
    title('Span Hat Vector Tips')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    hchords = figure ;
    hold on
    for j = 1:length(currt)
    plot3([spanHatR_xy(1,j) spanHatR_xy(1,j)+chordHatR_xy(1,j)],...
        [spanHatR_xy(2,j) spanHatR_xy(2,j)+chordHatR_xy(2,j)],...
        [spanHatR_xy(3,j) spanHatR_xy(3,j)+chordHatR_xy(3,j)],'Color','r')
    plot3([spanHatL_xy(1,j) spanHatL_xy(1,j)+chordHatL_xy(1,j)],...
        [spanHatL_xy(2,j) spanHatL_xy(2,j)+chordHatL_xy(2,j)],...
        [spanHatL_xy(3,j) spanHatL_xy(3,j)+chordHatL_xy(3,j)],'Color','b')
    end
    axis equal
    title('Chord Hat Vector Tips')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
end

%% Find vectors that give wing tip positions

span = .0025 ; %meters
AHat = repmat([0, 1/sqrt(2), 1/sqrt(2)],Nframes*repNum,1) ;
r_hinge = [0 .0002 0] ; % meters %[0 .0002 .0006]
bodyCM = zeros(Nframes*repNum,3) - (1/50e-6)*repmat(r_hinge,Nframes*repNum,1) ;

wingTipR = span*spanHatR_xy' ;
wingTipL = span*spanHatL_xy' ;


%% Save to data structure
clear('data')
data = struct ;
data.params.startTrackingTime = 0 ;
data.params.endTrackingTime = Nframes*repNum - 1 ;
data.params.fps = (currt(2)- currt(1))^(-1) ;
data.t = currt ;

data.backFlipTimesL = [] ;
data.backFlipTimesR = [] ;
data.fwdFlipTimesL = [] ;
data.fwdFlipTimesR = [] ;

for l = 1:repNum
    data.backFlipTimesL = [data.backFlipTimesL currt(Nframes*(l-1)+1)] ;
    data.backFlipTimesR = [data.backFlipTimesR currt(Nframes*(l-1)+1)] ;
    
    data.fwdFlipTimesL = [data.fwdFlipTimesL currt(Nframes*(l-1)+1+Nframes/2)] ;
    data.fwdFlipTimesR = [data.fwdFlipTimesL currt(Nframes*(l-1)+1+Nframes/2)] ;
end
data.rightWingTips = wingTipR/(50e-6) ; %convert to voxel units to match with 'data
data.leftWingTips = wingTipL/(50e-6) ;

data.rightSpanHats = spanHatR_xy' ;
data.leftSpanHats = spanHatL_xy' ;

data.rightChordHats = chordHatR_xy' ; 
data.leftChordHats = chordHatL_xy' ;

data.AHat = AHat ;
data.bodyCM = bodyCM ;

%% Find torques

run quasiSteadyTorqueSam
[htorque, havg, timeAvgTorque] = phiTorquePlot(data) ;

if saveFlag 
    cd(savePath)
    
    try
        mkdir(saveStr)
        cd(saveStr)
    catch
        cd(saveStr)
    end
    
    savefig(hkinematics,'kinematics')
    savefig(hspans, 'spans')
    savefig(hchords,'chords')
    savefig(htorque,'instantaneousTorques')
    savefig(havg,'avgTorques')
    save timeAvgTorque timeAvgTorque
end
    