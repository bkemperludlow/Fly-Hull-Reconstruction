function [t_latency, t_corrected] = ...
    findCorrectionTimes(data, bodyPitch, phiFrontAvg, t, pulseLength, pitchType)

%t_latency is time from 0 to when correction begins.
%t_corrected is time from 0 for fly to recover to within 10% of its
%original pitch angle

%pulseLength should be in ms, t in s
defineConstantsScript ;


zeropoint = find(t == 0) ;
[~, pulseEndInd] = min(abs(t - pulseLength/1000)) ; 
baseBodyPitch = bodyPitch(zeropoint) ;
bodyPitch_afterPert = bodyPitch(zeropoint:end) ;

phiR = -data.anglesBodyFrame(:,PHIR) ; 
phiL = data.anglesBodyFrame(:,PHIL) ; 

fwdFlipTimesR = data.fwdFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
%backFlipTimesR = data.backFlipTimesR ;
%backFlipTimesL = data.backFlipTimesL ;
phiFrontRs = zeros(size(fwdFlipTimesR)) ; 
phiFrontLs = zeros(size(fwdFlipTimesL)) ; 

for q = 1:length(fwdFlipTimesR)
    phiFrontRs(q) = phiR(find(t == fwdFlipTimesR(q))) ;
end
for s = 1:length(fwdFlipTimesL)
    phiFrontLs(s) = phiL(find(t == fwdFlipTimesL(s))) ;
end
phiDiffR = phiFrontRs - phiFrontAvg ; 
phiDiffL = phiFrontLs - phiFrontAvg ; 

if strcmp(pitchType,'up')
    [maxbodyPitch, maxind] = max(bodyPitch_afterPert) ;
    
    i1 = find(bodyPitch_afterPert <= .1*(maxbodyPitch-baseBodyPitch)+baseBodyPitch) ;
    i2 = find(i1 > maxind ) ;
    i3 = min(i1(i2)) ;
    
    t_corrected = t(i3+zeropoint-1) ; %seconds

    j1R = find(phiDiffR > 4) ;
    j1L = find(phiDiffL > 4) ; 
    
    if isempty(j1R) || isempty(j1L)
        t_latency = NaN ;
    else
        j2R = find(fwdFlipTimesR > pulseLength/1000) ;
        j2L = find(fwdFlipTimesL > pulseLength/1000) ;
        
        j3R = find(ismember(j2R,j1R) == 1, 1, 'first') ;
        j3L = find(ismember(j2L,j1L) == 1, 1, 'first') ;
        
        j4R = find(t == fwdFlipTimesR(j2R(j3R(1)))) ;
        j4L = find(t == fwdFlipTimesL(j2L(j3L(1)))) ;
        
        t_latency = .5*(t(j4R) + t(j4L)) ;
    end
elseif strcmp(pitchType,'down')
    [minbodyPitch, minind] = min(bodyPitch(1:800,1)) ;
    
    i1 = find(bodyPitch >= .1*(minbodyPitch-baseBodyPitch)+baseBodyPitch) ;
    i2 = find(i1 > minind ) ;
    i3 = min(i1(i2)) ;
    
    t_corrected = t(i3) ; %seconds
    
    j1R = find(phiDiffR < 4) ;
    j1L = find(phiDiffL < 4) ;
    
    if isempty(j1R) || isempty(j1L)
        t_latency = NaN ;
    else
        
        j2R = find(fwdFlipTimesR > pulseLength/1000) ;
        j2L = find(fwdFlipTimesL > pulseLength/1000) ;
        
        j3R = find(ismember(j2R,j1R) == 1, 1, 'first') ;
        j3L = find(ismember(j2L,j1L) == 1, 1, 'first') ;
        
        j4R = find(t == fwdFlipTimesR(j2R(j3R(1))) ) ;
        j4L = find(t == fwdFlipTimesL(j2L(j3L(1))) ) ;
        
        t_latency = .5*(t(j4R) + t(j4L)) ;
    end
else
    disp('check pitchType')
    return ;
end