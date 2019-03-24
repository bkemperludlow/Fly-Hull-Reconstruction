%function [h_main, phiFront, fwdFlipTimes] = plotPitchAndPhiFront(ExprNum, MovNum) 

ExprNum = 7 ;
MovNum = 8 ;

if MovNum < 10
    zstr = '00' ;
elseif MovNum < 100 ;
    zstr = '0' ;
else
    zstr = '' ;
end

movieStr = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum)] ;
dataPath = ['G:\Janelia Flies\kir2.1 flies\Analysis\No Perturbation\' movieStr] ;
cd(dataPath)
fileName = [movieStr '_test.mat'] ;
load(fileName) ;

defineConstantsScript

fwdFlipTimesR = data.fwdFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;

t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;
%{
fwdFlipIndR = zeros(size(fwdFlipTimesR)) ;
fwdFlipIndL = zeros(size(fwdFlipTimesL)) ;

for i = 1:length(fwdFlipTimesR)
    fwdFlipIndR(i) = find(t == fwdFlipTimesR(i)) ;
end
for j = 1:length(fwdFlipTimesL)
    fwdFlipIndL(j) = find(t == fwdFlipTimesL(j)) ;
end
%}
phiR = -data.anglesBodyFrame(:,PHIR) ;
phiL = data.anglesBodyFrame(:,PHIL) ;
bodyPitch = data.anglesLabFrame(:,BETA) ;

phiEstErr = 2 ;
[sp_phiR, ~, ~] = mySplineSmooth(t,phiR,phiEstErr) ;
[sp_phiL, ~, ~] = mySplineSmooth(t,phiL,phiEstErr) ;

pitchEstErr = .75 ;
[sp_pitch, ~,~] = mySplineSmooth(t,bodyPitch,pitchEstErr) ;
pitch_smoothed = fnval(sp_pitch, t) ;

if (0)
    figure ; hold on
    plot(t*1000, phiR,'ro')
    plot(t*1000, phiL,'bo')
    plot(t*1000, fnval(sp_phiR,t) ,'r-')
    plot(t*1000, fnval(sp_phiL,t),'b-')
    xlabel('Time [ms]')
    ylabel('Stroke Angle, \Phi , [deg]')
    
    axis tight
    
    figure ; hold on
    plot(t*1000, bodyPitch, 'k.')
    plot(t*1000, pitch_smoothed, 'r-')
    xlabel('Time [ms]')
    ylabel('Body Pitch, \theta_b , [deg]')
    axis tight
end

phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;

if length(phiFrontR) == length(phiFrontL)
    phiFront = (phiFrontR + phiFrontL) / 2 ;
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL ) /2 ;
elseif length(fwdFlipTimesR) < length(fwdFlipTimesL)
    idx = zeros(length(fwdFlipTimesR),1) ;
    for q = 1:length(fwdFlipTimesR)
        [~,minInd] = min(abs(fwdFlipTimesL - fwdFlipTimesR(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL(idx)) / 2 ;
    phiFront = (phiFrontR + phiFrontL(idx)) / 2 ;
elseif length(fwdFlipTimesL) < length(fwdFlipTimesR)
    idx = zeros(length(fwdFlipTimesL),1) ;
    for q = 1:length(fwdFlipTimesL)
        [~,minInd] = min(abs(fwdFlipTimesR - fwdFlipTimesL(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesL + fwdFlipTimesR(idx)) / 2 ;
    phiFront = (phiFrontL + phiFrontR(idx)) / 2 ;
end

bodyPitch_meanSubtracted = pitch_smoothed - mean(pitch_smoothed) ;
phiFront_meanSubtracted = phiFront - mean(phiFront) ;

bodyPitch_MS_downSampled = fnval(sp_pitch, fwdFlipTimes) - mean(pitch_smoothed) ;

figure 
hold on
plot(fwdFlipTimes*1000, bodyPitch_MS_downSampled, 'r-')
plot(fwdFlipTimes*1000, -phiFront_meanSubtracted, 'b-')
title(['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum)])



[cor, lag] = xcorr(-phiFront_meanSubtracted,bodyPitch_MS_downSampled) ;
figure ;
plot(lag, cor)
title(['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum)])



