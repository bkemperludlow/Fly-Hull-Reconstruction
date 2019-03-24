function [fwdFlipTimesR, backFlipTimesR, fwdFlipTimesL, backFlipTimesL, anglesBodyFrame,...
    anglesLabFrame, data] = saveWingFlipsAndAngles(exprNum,movNum,PitchType)
% Not all of the 'data' files are consistent, because the code is often
% updated, so this function calculates and saves the flip times and
%(unsmoothed) body angles. This way I can just throw it in an 'if'
%statement if my data structure doesn't have what I want
%
%Note that this is just the first half of quick_and_dirty.m, minus the
%plots
%
%PitchType is the only ambiguous variable. It should have value either 'up'
%or 'down' (string)
%---------------------------------------------------------------------------------
defineConstantsScript

if exprNum == 7
    pulseLengthMS = 5.8 ;
else 
    pulseLengthMS = 8 ;
end
pulseStartMS = 0 ;

if (movNum<10)
    zstr = '00' ;
elseif (movNum<100)
    zstr = '0' ;
else
    zstr = '';
end

if strcmp(PitchType,'up')
    datapath = ['F:\luca\Analysis\pitch up\Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '\' ] ;
elseif strcmp(PitchType,'down')
    datapath = ['F:\luca\Analysis\pitch down\Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '\' ] ;
else
    disp('check PitchType')
    return ;
end

datafilename = [ datapath ...
    'Expr' num2str(exprNum) 'mov' zstr num2str(movNum) '_Data_manually_corrected.mat' ] ; %

load(datafilename) ;

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
rollEstErr = .5;
[anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody,sp_rho, smoothed_rho, rho_t, rho_samp] ...
    = calcAngles_quick_and_dirty_mk2(data, rollEstErr, false) ;

% spline smooth phiL and phiR
ignoreIndR = unique([find(isnan(anglesBodyFrame(:, PHIR))==1)' ignoreFrames]) ;% [127 311 312 425 908 909] ; % expr7mov9
ignoreIndL = unique([find(isnan(anglesBodyFrame(:, PHIL))==1)'  ignoreFrames])  ; % [425 536 537 869 908 1061] ; % expr7mov9

phiR = -anglesBodyFrame(:, PHIR) ;
if (~isempty(ignoreIndR))
    phiR(ignoreIndR) = NaN ;
end

phiL = +anglesBodyFrame(:, PHIL) ;
if (~isempty(ignoreIndL))
    phiL(ignoreIndL) = NaN ;
end

[fwdFlipTimesR, backFlipTimesR, fwdFlipIndR, backFlipIndR, fwdFlipPhiR, backFlipPhiR, badIndicesR] = findWingFlipTimes_mk3 (t, phiR, false);
[fwdFlipTimesL, backFlipTimesL, fwdFlipIndL, backFlipIndL, fwdFlipPhiL, backFlipPhiL, badIndicesL] = findWingFlipTimes_mk3 (t, phiL, false);

phiR(badIndicesR) = NaN ;
phiL(badIndicesL) = NaN ;

data.fwdFlipTimesR = fwdFlipTimesR ;
data.backFlipTimesR = backFlipTimesR ;
data.fwdFlipTimesL = fwdFlipTimesL ;
data.backFlipTimesL = backFlipTimesL ;
data.params.pulseLengthMS = pulseLengthMS ; 
data.anglesBodyFrame = anglesBodyFrame ;
data.anglesLabFrame = anglesLabFrame ;

save(datafilename,'data')