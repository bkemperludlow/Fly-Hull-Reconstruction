function [fwdStrokes, backStrokes] = analyzeHalfStrokes(data, plotFlag) 

%plotFlag = false ;

%% load data
fwdFlipTimesR = data.fwdFlipTimesR ;
backFlipTimesR = data.backFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
backFlipTimesL = data.backFlipTimesL ;

if isfield(data,'manualCorrRangeMS')
    t_range = data.manualCorrRangeMS / 1000 ;
else
    t_range = [data.params.startTrackingTime , data.params.endTrackingTime ] / 8000 ;    
end
tms = (data.params.startTrackingTime : data.params.endTrackingTime) / 8 ;

%% check kinematics
if plotFlag
    defineConstantsScript
    patchColor = [1 1 1 ] * 0.8;
    phiR = -data.anglesBodyFrame(:,PHIR) ;
    phiL = data.anglesBodyFrame(:,PHIL) ; 
    figure ;
    hold on
    plot(tms, phiR,'ro-')
    plot(tms, phiL,'bo-')
    plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
    set(gca,'xlim',t_range*1000)
    set(gca,'ylim',[0 200])
end

%% sort out forward flip times
if length(fwdFlipTimesR) == length(fwdFlipTimesL)
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL ) /2 ;
elseif length(fwdFlipTimesR) < length(fwdFlipTimesL)
    idx = zeros(length(fwdFlipTimesR),1) ;
    for q = 1:length(fwdFlipTimesR)
        [~,minInd] = min(abs(fwdFlipTimesL - fwdFlipTimesR(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL(idx)) / 2 ;
elseif length(fwdFlipTimesL) < length(fwdFlipTimesR)
    idx = zeros(length(fwdFlipTimesL),1) ;
    for q = 1:length(fwdFlipTimesL)
        [~,minInd] = min(abs(fwdFlipTimesR - fwdFlipTimesL(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesL + fwdFlipTimesR(idx)) / 2 ;
end

%% sort out backward flip times
if length(backFlipTimesR) == length(backFlipTimesL)
    backFlipTimes = (backFlipTimesR + backFlipTimesL ) /2 ;
elseif length(backFlipTimesR) < length(backFlipTimesL)
    idx = zeros(length(backFlipTimesR),1) ;
    for q = 1:length(backFlipTimesR)
        [~,minInd] = min(abs(backFlipTimesL - backFlipTimesR(q))) ; 
        idx(q) = minInd ;
    end
    backFlipTimes = (backFlipTimesR + backFlipTimesL(idx)) / 2 ;
elseif length(backFlipTimesL) < length(backFlipTimesR)
    idx = zeros(length(backFlipTimesL),1) ;
    for q = 1:length(backFlipTimesL)
        [~,minInd] = min(abs(backFlipTimesR - backFlipTimesL(q))) ; 
        idx(q) = minInd ;
    end
    backFlipTimes = (backFlipTimesL + backFlipTimesR(idx)) / 2 ;
end

%% restrict analysis to regions in which kinematics are (hopefully) trustworthy
fwdFlipTimes = fwdFlipTimes(fwdFlipTimes > t_range(1) & fwdFlipTimes < t_range(2)) ; 
backFlipTimes = backFlipTimes(backFlipTimes > t_range(1) & backFlipTimes < t_range(2)) ; 

%% find duration of forward strokes (fwdFlip(t+1) - backFlip(t))

N_fwdFlips = length(fwdFlipTimes) ;
fwdStrokes = nan(N_fwdFlips,1) ;
for i = 1:N_fwdFlips
    fwdFlip = fwdFlipTimes(i) ;
    ind_temp = find(backFlipTimes < fwdFlip, 1, 'last') ;
    if ~isempty(ind_temp)
        fwdStrokes(i) = fwdFlip - backFlipTimes(ind_temp) ; 
    end
end
fwdStrokes = fwdStrokes(~isnan(fwdStrokes)) ;

%% find duration of backward strokes (backFlip(t+1) - fwdFlip(t))

N_backFlips = length(backFlipTimes) ;
backStrokes = nan(N_backFlips,1) ;
for i = 1:N_backFlips
    backFlip = backFlipTimes(i) ;
    ind_temp = find(fwdFlipTimes < backFlip, 1, 'last') ;
    if ~isempty(ind_temp)
        backStrokes(i) = backFlip - fwdFlipTimes(ind_temp) ; 
    end
end
backStrokes = backStrokes(~isnan(backStrokes)) ;


end