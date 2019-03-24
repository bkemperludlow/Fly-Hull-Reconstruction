function [phiFront, phiBack, phiAmp, phiMid, rawPhiFront, rawPhiBack, ...
    rawPhiAmp, rawPhiMid] = strokeAngleChange(data,polyCoeffs,t)

%To get back the avg ones, just take e.g. rawPhiFront - phiFront

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

phiR = (-1)*data.anglesBodyFrame(:,PHIR) ;
phiL = data.anglesBodyFrame(:,PHIL) ;

%prePertInd = find((t > manualCorrRangeS(1)) & (t < 0)) ;
%avgDeltaPhi =  ;

prePertFwdFlipIndR = find((fwdFlipTimesR > manualCorrRangeS(1)) & (fwdFlipTimesR < 0)) ;
prePertFwdFlipIndL = find((fwdFlipTimesL > manualCorrRangeS(1)) & (fwdFlipTimesL < 0)) ;
prePertFwdIndR = zeros(size(prePertFwdFlipIndR)) ;
prePertFwdIndL = zeros(size(prePertFwdFlipIndL)) ;
for j = 1:length(prePertFwdFlipIndR)
    prePertFwdIndR(j) = find(t == fwdFlipTimesR(prePertFwdFlipIndR(j))) ;
end
for k = 1:length(prePertFwdFlipIndL)
    prePertFwdIndL(k) = find(t == fwdFlipTimesL(prePertFwdFlipIndL(k))) ;
end
prePertPhiFrontR = phiR(prePertFwdIndR) ;
prePertPhiFrontL = phiL(prePertFwdIndL) ;
avgPhiFront = mean([mean(prePertPhiFrontR) ; mean(prePertPhiFrontL)]) ;

prePertBackFlipIndR = find((backFlipTimesR > manualCorrRangeS(1)) & (backFlipTimesR < 0)) ;
prePertBackFlipIndL = find((backFlipTimesL > manualCorrRangeS(1)) & (backFlipTimesL < 0)) ;
prePertBackIndR = zeros(size(prePertBackFlipIndR)) ;
prePertBackIndL = zeros(size(prePertBackFlipIndL)) ;
for m = 1:length(prePertBackFlipIndR)
    prePertBackIndR(m) = find(t == backFlipTimesR(prePertBackFlipIndR(m))) ;
end
for n = 1:length(prePertBackFlipIndL)
    prePertBackIndL(n) = find(t == backFlipTimesL(prePertBackFlipIndL(n))) ; 
end
prePertPhiBackR = phiR(prePertBackIndR) ;
prePertPhiBackL = phiL(prePertBackIndL) ;
avgPhiBack = mean([mean(prePertPhiBackR) ; mean(prePertPhiBackL)]) ;

avgPhiMid = (avgPhiBack+avgPhiFront)/2 ;
avgPhiAmp = avgPhiBack - avgPhiFront ;

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

rawPhiFront = mean([phiR(fwdMinIndR) phiL(fwdMinIndL)]) ;
rawPhiBack =mean([phiR(backMinIndR) phiL(backMinIndL)]) ;
rawPhiMid = (rawPhiBack + rawPhiFront)/2 ;
rawPhiAmp = rawPhiBack - rawPhiFront ;


phiFront = rawPhiFront - avgPhiFront ;
phiBack = rawPhiBack - avgPhiBack ;
phiAmp = rawPhiAmp - avgPhiAmp ;
phiMid = rawPhiMid - avgPhiMid ; 



