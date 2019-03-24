function [etaFwd, etaBack, rawEtaFwd, rawEtaBack] = wingpitchAngleChange(data,polyCoeffs,t,plotFlag)

%eta is wing pitch angle. etaMax and etaMin refer to deviations from
%'normal', rawEtaMax and rawEtaMin are just the values of eta averaged over
%a halfstroke
%
%Fwd should be minimum, and back should be maximum
if ~exist('plotFlag','var')
    plotFlag = 0 ;
end

defineConstantsScript

backFlipTimesR = data.backFlipTimesR ;
fwdFlipTimesR = data.fwdFlipTimesR ;
backFlipTimesL = data.backFlipTimesL ;
fwdFlipTimesL = data.fwdFlipTimesL ;
if (isfield(data,'correctionTime'))
    manualCorrRangeS =data.correctionTime/1000 ;
elseif(isfield(data,'manualCorrRangeMS'))
    manualCorrRangeS =data.manualCorrRangeMS/1000 ;
else
    manualCorrRangeS = [-10 30] ; 
end

etaR = (180/pi)*unwrap((pi/180)*data.anglesBodyFrame(:,ETAR)) ;
etaL = (180/pi)*unwrap((pi/180)*data.anglesBodyFrame(:,ETAL)) ;

%correct for errors in eta
tooHighEtaR = find(etaR > 180) ; 
etaR(tooHighEtaR) = etaR(tooHighEtaR) - 360 ; 
tooHighEtaL = find(etaL > 180) ; 
etaL(tooHighEtaL) = etaR(tooHighEtaL) - 360 ;
tooLowEtaR = find(etaR < 0) ;
etaR(tooLowEtaR) = etaR(tooLowEtaR) + 360 ; 
tooLowEtaL = find(etaL < 0) ;
etaL(tooLowEtaL) = etaL(tooLowEtaL) + 360 ; 

%Find pre-perturbation indices
prePertFwdFlipIndR = find((fwdFlipTimesR > manualCorrRangeS(1)) & (fwdFlipTimesR < 0)) ;
prePertFwdFlipIndL = find((fwdFlipTimesL > manualCorrRangeS(1)) & (fwdFlipTimesL < 0)) ;
prePertFwdIndR = ones(length(prePertFwdFlipIndR),2) ;
prePertFwdIndL = ones(length(prePertFwdFlipIndL),2) ;
for j = 1:length(prePertFwdFlipIndR)
    prePertFwdIndR(j,1) = find(t == fwdFlipTimesR(prePertFwdFlipIndR(j))) ;
end
for k = 1:length(prePertFwdFlipIndL)
    prePertFwdIndL(k,1) = find(t == fwdFlipTimesL(prePertFwdFlipIndL(k))) ;
end

prePertBackFlipIndR = find((backFlipTimesR > manualCorrRangeS(1)) & (backFlipTimesR < 0)) ;
prePertBackFlipIndL = find((backFlipTimesL > manualCorrRangeS(1)) & (backFlipTimesL < 0)) ;
prePertBackIndR = zeros(length(prePertBackFlipIndR),2) ;
prePertBackIndL = zeros(length(prePertBackFlipIndL),2) ;
for m = 1:length(prePertBackFlipIndR)
    prePertBackIndR(m,1) = find(t == backFlipTimesR(prePertBackFlipIndR(m))) ;
end
for n = 1:length(prePertBackFlipIndL)
    prePertBackIndL(n,1) = find(t == backFlipTimesL(prePertBackFlipIndL(n))) ; 
end

prePertIndR = sortrows([prePertBackIndR;prePertFwdIndR],1) ;
prePertIndL = sortrows([prePertBackIndL;prePertFwdIndL],1) ;

EtaRAvgs = zeros(size(prePertIndR,1)-1, 2) ;
EtaLAvgs = zeros(size(prePertIndL,1)-1, 2) ;

for q = 1:size(EtaRAvgs,1)
    EtaRAvgs(q,1) = mean(etaR(prePertIndR(q,1):prePertIndR(q+1,1))) ; 
    EtaRAvgs(q,2) = prePertIndR(q+1,2) ; 
end
for p = 1:size(EtaLAvgs,1)
    EtaLAvgs(p,1) = mean(etaL(prePertIndL(p,1):prePertIndL(p+1,1))) ; 
    EtaLAvgs(p,2) = prePertIndL(p+1,2) ; 
end

EtaRFwdInd = find(EtaRAvgs(:,2) == 1) ;
EtaRBackInd = find(EtaRAvgs(:,2) == 0) ;
etaRFwdAvg = mean(EtaRAvgs(EtaRFwdInd,1)) ;
etaRBackAvg = mean(EtaRAvgs(EtaRBackInd,1)) ;

EtaLFwdInd = find(EtaLAvgs(:,2) == 1) ;
EtaLBackInd = find(EtaLAvgs(:,2) == 0) ;
etaLFwdAvg = mean(EtaLAvgs(EtaLFwdInd,1)) ;
etaLBackAvg = mean(EtaLAvgs(EtaLBackInd,1)) ;

etaFwdAvg = (etaLFwdAvg + etaRFwdAvg)/2 ;
etaBackAvg = (etaLBackAvg + etaRBackAvg)/2 ;

peakAccel = -polyCoeffs(2)/(2*polyCoeffs(1)) ;
backFlipDiffR = abs(backFlipTimesR - peakAccel) ;
fwdFlipDiffR = abs(fwdFlipTimesR - peakAccel) ;
backFlipDiffL = abs(backFlipTimesL - peakAccel) ;
fwdFlipDiffL = abs(fwdFlipTimesL - peakAccel) ;
[~, BMRInd] = min(backFlipDiffR) ;
[~, BMLInd] = min(backFlipDiffL) ;
[~, FMRInd] = min(fwdFlipDiffR) ;
[~, FMLInd] = min(fwdFlipDiffL) ;

backMinIndR = find(t == backFlipTimesR(BMRInd)) ; 
backMinIndL = find(t == backFlipTimesL(BMLInd)) ; 
fwdMinIndR = find(t == fwdFlipTimesR(FMRInd)) ; 
fwdMinIndL = find(t == fwdFlipTimesL(FMLInd)) ;

if backMinIndR < fwdMinIndR
    rawEtaFwdR = mean(etaR(backMinIndR:fwdMinIndR)) ; 
    rawEtaFwdL = mean(etaL(backMinIndL:fwdMinIndL)) ;
    rawEtaFwd = (rawEtaFwdR + rawEtaFwdL)/2 ; 
    
    fwdFlipBeforeIndR = find(t == fwdFlipTimesR(FMRInd-1)) ;
    fwdFlipBeforeIndL = find(t == fwdFlipTimesL(FMLInd-1)) ;
    
    rawEtaBackR = mean(etaR(fwdFlipBeforeIndR:backMinIndR)) ; 
    rawEtaBackL = mean(etaL(fwdFlipBeforeIndL:backMinIndL)) ;
    rawEtaBack = (rawEtaBackR + rawEtaBackL)/2 ;
elseif fwdMinIndR < backMinIndR 
    rawEtaBackR = mean(etaR(fwdMinIndR:backMinIndR)) ; 
    rawEtaBackL = mean(etaL(fwdMinIndL:backMinIndL)) ;
    rawEtaBack = (rawEtaBackR + rawEtaBackL)/2 ; 
    
    backFlipBeforeIndR = find(t == backFlipTimesR(BMRInd-1)) ;
    backFlipBeforeIndL = find(t == backFlipTimesL(BMLInd-1)) ;
    
    rawEtaFwdR = mean(etaR(backFlipBeforeIndR:fwdMinIndR)) ; 
    rawEtaFwdL = mean(etaL(backFlipBeforeIndL:fwdMinIndL)) ;
    rawEtaFwd = (rawEtaFwdR + rawEtaFwdL)/2 ;
else 
    disp('check something')
    return ;
end

etaFwd = rawEtaFwd - etaFwdAvg ;
etaBack = rawEtaBack - etaBackAvg ; 

if plotFlag
    heta = figure ;
    phiR = data.anglesBodyFrame(:,PHIR) ;
    phiL = data.anglesBodyFrame(:,PHIL) ;
    
    s1 = subplot(2,2,1) ;
        plot(phiL(prePertIndL),etaL(prePertIndL),'b','LineWidth',2)
        xlabel('\phi_L [deg]')
        ylabel('\eta_L [deg]')
        title('Left Wing, Pre-Perturbation')
    s2 = subplot(2,2,2) ;
        plot(-phiR(prePertIndR),etaR(prePertIndR),'r','LineWidth',2)
        xlabel('\phi_R [deg]')
        ylabel('\eta_R [deg]')
        title('Right Wing, Pre-Perturbation')
    s3 = subplot(2,2,3) ;
        %if 
        plot(phiL(prePertIndR),etaL(prePertIndR),'b','LineWidth',2)
        xlabel('\phi_L [deg]')
        ylabel('\eta_L [deg]')
        title('Left Wing, Post-Perturbation')
    s4 = subplot(2,2,4) ;
        plot(-phiR(prePertIndR),etaR(prePertIndR),'r','LineWidth',2)
        xlabel('\phi_R [deg]')
        ylabel('\eta_R [deg]')
        title('Right Wing, Post-Perturbation')
end


    
    
    
    
