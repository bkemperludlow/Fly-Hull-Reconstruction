ExprNum =  [ 7 7 7 7 7 7 13] ;%[7 7 7 7 14] ; 
    %[7 7 7] ; 
MovNum = {'007' '008' '009' '015' '054' '062' '007'};  %{'006' '040' '065' '067' '023'} ;%
    %{'040' '065' '067'} ;
runSize = length(ExprNum) ;
defineConstantsScript 
velocityThreshold = 100 ; %deg/s

Mov = zeros(runSize,1) ;
Expr = ExprNum' ;
T_back = zeros(runSize,1) ;
T_front = zeros(runSize,1) ;
T_mid_front = zeros(runSize,1) ;
T_mid_back = zeros(runSize,1) ;

for i = 1:runSize
    datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
            num2str(ExprNum(i)), '_mov_',MovNum{i}) ;
    datafilename = strcat(datapath,'\Expr',num2str(ExprNum(i)), ...
        'mov',MovNum{i}, '_Data_manually_corrected.mat') ;
    
    load(datafilename) ;
    Mov(i) = str2num(MovNum{i}) ;
    
    fwdFlipTimesR = data.fwdFlipTimesR ; 
    backFlipTimesR = data.backFlipTimesR ; 
    %fwdFlipTimesL = data.fwdFlipTimesL ; 
    %backFlipTimesL = data.backFlipTimesL ; 
    
    startTime = data.params.startTrackingTime ;
    endTime = data.params.endTrackingTime ;
    fps = data.params.fps ;
    t = (startTime:endTime)/fps ; %seconds
    
    bodyPitch = data.anglesLabFrame(:,BETA) ;
    EstErr = 1 ;
    [sp_pitch, ~,~] = mySplineSmooth(t,bodyPitch,EstErr) ; 
    pitchVelocity = fnder(sp_pitch,1) ;
    
    backInd = find(backFlipTimesR > 0, 1, 'first') ;
    while abs(fnval(pitchVelocity, backFlipTimesR(backInd))) < velocityThreshold
        backInd = backInd + 1 ;
    end
    T_back(i) = backFlipTimesR(backInd) ;
    
    fwdInd = find(fwdFlipTimesR > 0, 1, 'first') ;
    while abs(fnval(pitchVelocity, fwdFlipTimesR(fwdInd))) < velocityThreshold
        fwdInd = fwdInd + 1 ;
    end
    T_front(i) = fwdFlipTimesR(fwdInd) ;
    
    if fwdFlipTimesR(fwdInd) < backFlipTimesR(backInd)
        temp1 = (fwdFlipTimesR(fwdInd) + backFlipTimesR(backInd-1))/2 ;
        temp2 = (fwdFlipTimesR(fwdInd+1) + backFlipTimesR(backInd))/2 ;
        if (fnval(pitchVelocity, temp1) > velocityThreshold) && (temp1 > 0)
            T_mid_front(i) = temp1 ;
        else
            T_mid_front(i) = temp2 ;
        end
        T_mid_back(i) = (fwdFlipTimesR(fwdInd) + backFlipTimesR(backInd))/2 ;
    elseif backFlipTimesR(backInd) < fwdFlipTimesR(fwdInd) 
        temp1 = (fwdFlipTimesR(fwdInd-1) + backFlipTimesR(backInd))/2 ;
        temp2 = (fwdFlipTimesR(fwdInd) + backFlipTimesR(backInd+1))/2 ;
        if (fnval(pitchVelocity, temp1) > velocityThreshold) && (temp1 > 0)
            T_mid_back(i) = temp1 ;
        else
            T_mid_back(i) = temp2 ;
        end
        T_mid_front(i) = (fwdFlipTimesR(fwdInd) + backFlipTimesR(backInd))/2 ;
    else
        disp('ERROR')
        return ;
    end
        
end

PhaseData = table(Expr, Mov, T_back, T_front, T_mid_back, T_mid_front) ; 
