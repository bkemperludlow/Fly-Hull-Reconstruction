%% Preliminaries
cd('F:\luca\Analysis\pitch up\Expr_7_mov_008')
load('Expr7mov008_Data_manually_corrected.mat')
defineConstantsScript

phiCut = false ;
etaCut = false ;
thetaCut = false ;

EstPhiErr = 1 ;
EstEtaErr = 1 ; %2.5
EstThetaErr = 1 ; %2

checkPlotFlag1 = false ;

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

phiR = -data.anglesBodyFrame(:,PHIR) ;
etaR = (180/pi)*unwrap((pi/180)*data.anglesBodyFrame(:,ETAR)) ;
thetaR = data.anglesBodyFrame(:,THETAR) ;
phiR(ignoreFrames) = NaN ; etaR(ignoreFrames) = NaN ; thetaR(ignoreFrames) = NaN ;

phiL = data.anglesBodyFrame(:,PHIL) ;
etaL = (180/pi)*unwrap((pi/180)*data.anglesBodyFrame(:,ETAL)) ;
thetaL = data.anglesBodyFrame(:,THETAL) ;
phiL(ignoreFrames) = NaN ; etaL(ignoreFrames) = NaN ; thetaL(ignoreFrames) = NaN ;

bodyCM = data.bodyCM ;
thetaB = data.anglesLabFrame(:,THETAB) ;
rhoB = data.anglesLabFrame(:,RHO) ;
psiB = data.anglesLabFrame(:,PSIB) ;

%% Find times interval to cut and paste (cut [cutpoint:zeropoint] and paste afterwards)
zeropoint = find(t==0) ;
prePertFlipInds = find(backFlipTimesR < 0) ;
avgPeriod = mean(diff(backFlipTimesR(prePertFlipInds))) ;
[~,cutpoint] = min(abs(t + avgPeriod)) ;

if phiCut
    phiR_new = [phiR(1:zeropoint) ; repmat(phiR(cutpoint:zeropoint),20,1)] ;
    phiL_new = [phiL(1:zeropoint) ; repmat(phiL(cutpoint:zeropoint),20,1)] ;
    Ntemp = length(t) ;
    phiR_new = phiR_new(1:Ntemp) ;
    phiL_new = phiL_new(1:Ntemp) ;
else 
    phiR_new = phiR ;
    phiL_new = phiL ;
end

if etaCut
    etaR_new = [etaR(1:zeropoint) ; repmat(etaR(cutpoint:zeropoint),20,1)] ;
    etaL_new = [etaL(1:zeropoint) ; repmat(etaL(cutpoint:zeropoint),20,1)] ;
    Ntemp = length(t) ;
    etaR_new = etaR_new(1:Ntemp) ;
    etaL_new = etaL_new(1:Ntemp) ;
else 
    etaR_new = etaR ;
    etaL_new = etaL ;
end

if thetaCut
    thetaR_new = [thetaR(1:zeropoint) ; repmat(thetaR(cutpoint:zeropoint),20,1)] ;
    thetaL_new = [thetaL(1:zeropoint) ; repmat(thetaL(cutpoint:zeropoint),20,1)] ;
    Ntemp = length(t) ;
    thetaR_new = thetaR_new(1:Ntemp) ;
    thetaL_new = thetaL_new(1:Ntemp) ;
else 
    thetaR_new = thetaR ;
    thetaL_new = thetaL ;
end

%% Smooth new angle data
[sp_phiR, ~, ~] = mySplineSmooth(t,phiR_new,EstPhiErr) ;
[sp_phiL, ~, ~] = mySplineSmooth(t,phiL_new,EstPhiErr) ;

[sp_etaR, ~, ~] = mySplineSmooth(t,etaR_new,EstEtaErr) ;
[sp_etaL, ~, ~] = mySplineSmooth(t,etaL_new,EstEtaErr) ;

[sp_thetaR, ~, ~] = mySplineSmooth(t,thetaR_new,EstThetaErr) ;
[sp_thetaL, ~, ~] = mySplineSmooth(t,thetaL_new,EstThetaErr) ;

phiR_smooth = fnval(sp_phiR,t)' ;
phiL_smooth = fnval(sp_phiL,t)' ;
etaR_smooth = fnval(sp_etaR,t)' ;
etaL_smooth = fnval(sp_etaL,t)' ;
thetaR_smooth = fnval(sp_thetaR,t)' ;
thetaL_smooth = fnval(sp_thetaL,t)' ;

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

%% Define wing vectors based on Euler angles

%these are the vector hats for the wings assuming that the stroke plane is
%coincident with the xy plane and the body axis is pointing along the y
%direction
spanHatR_xy = [sind(phiR_smooth).*cosd(thetaR_smooth); cosd(phiR_smooth).*cosd(thetaR_smooth); sind(thetaR_smooth)] ;
spanHatL_xy = [-sind(phiL_smooth).*cosd(thetaL_smooth); cosd(phiL_smooth).*cosd(thetaL_smooth); sind(thetaL_smooth)] ;

chordHatR_xy = [-cosd(etaR_smooth).*cosd(phiR_smooth); cosd(etaR_smooth).*sind(phiR_smooth) ; sind(etaR_smooth) ] ;
chordHatL_xy = [cosd(etaL_smooth).*cosd(phiL_smooth); cosd(etaL_smooth).*sind(phiL_smooth) ; sind(etaL_smooth) ] ;

badIndR = find(spanHatR_xy(1,:) < 0) ; 
badIndL = find(spanHatL_xy(1,:) > 0) ;
badInd = unique([badIndR, badIndL]) ;

NanMat = repmat([NaN; NaN; NaN], 1, length(badInd)) ;
spanHatR_xy(:,badInd) = NanMat ;
spanHatL_xy(:,badInd) = NanMat ;
chordHatR_xy(:,badInd) = NanMat ;
chordHatL_xy(:,badInd) = NanMat ;

if checkPlotFlag1
    figure; 
    idx1 = 840 ;
    idx2 = 875 ;
    plot3(spanHatR_xy(1,idx1:idx2),spanHatR_xy(2,idx1:idx2),spanHatR_xy(3,idx1:idx2),'ko','MarkerFaceColor','r')
    hold on
    plot3(spanHatL_xy(1,idx1:idx2),spanHatL_xy(2,idx1:idx2),spanHatL_xy(3,idx1:idx2),'ko','MarkerFaceColor','b')
    axis equal
    
    figure ;
    hold on
    idx1 = 740 ;
    idx2 = 780 ;
    spacing = 1 ;
    for j = idx1:spacing:idx2
    plot3([spanHatR_xy(1,j) spanHatR_xy(1,j)+chordHatR_xy(1,j)],...
        [spanHatR_xy(2,j) spanHatR_xy(2,j)+chordHatR_xy(2,j)],...
        [spanHatR_xy(3,j) spanHatR_xy(3,j)+chordHatR_xy(3,j)],'Color','r')
    plot3([spanHatL_xy(1,j) spanHatL_xy(1,j)+chordHatL_xy(1,j)],...
        [spanHatL_xy(2,j) spanHatL_xy(2,j)+chordHatL_xy(2,j)],...
        [spanHatL_xy(3,j) spanHatL_xy(3,j)+chordHatL_xy(3,j)],'Color','b')
    end
    axis equal
    
end

%% Find vectors that give wing tip positions

span = .0025 ; %meters
AHat = repmat([0, 1/sqrt(2), 1/sqrt(2)],Nframes,1) ;
r_hinge = [0 .0002 .0006] ; % meters
bodyCM = zeros(Nframes,3) - (1/50e-6)*repmat(r_hinge,Nframes,1) ;

wingTipR = span*spanHatR_xy' ;
wingTipL = span*spanHatL_xy' ;


%% Save to data structure
data.rightWingTips = wingTipR/(50e-6) ; %convert to voxel units to match with 'data
data.leftWingTips = wingTipL/(50e-6) ;

data.rightSpanHats = spanHatR_xy' ;
data.leftSpanHats = spanHatL_xy' ;

data.rightChordHats = chordHatR_xy' ; 
data.leftChordHats = chordHatL_xy' ;

data.AHat = AHat ;
data.bodyCM = bodyCM ;



    