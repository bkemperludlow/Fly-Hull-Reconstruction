ExprNum =  [ 7 7 7 7 7 7 7 7 7 7 13 14] ;%[7 7 7 7 14] ;
MovNum = {'007' '008' '009' '015' '019' '023' '024' '046' '054' '062' '007' '082'};  %{'006' '040' '065' '067' '023'} ;%

runSize = length(ExprNum) ;

Mov = zeros(runSize,1) ;
Expr = ExprNum' ;
pulsePhaseArrayR = zeros(runSize,1) ;
pulsePhaseArrayL = zeros(runSize,1) ;

for i = 1:runSize
    datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
            num2str(ExprNum(i)), '_mov_',MovNum{i}) ;
    datafilename = strcat(datapath,'\Expr',num2str(ExprNum(i)), ...
        'mov',MovNum{i}, '_Data_manually_corrected.mat') ;
    
    load(datafilename) ;
    Mov(i) = str2num(MovNum{i}) ;
    
    fwdFlipTimesR = data.fwdFlipTimesR ; 
    backFlipTimesR = data.backFlipTimesR ; 
    fwdFlipTimesL = data.fwdFlipTimesL ; 
    backFlipTimesL = data.backFlipTimesL ; 
    
    startTime = data.params.startTrackingTime ;
    endTime = data.params.endTrackingTime ;
    fps = data.params.fps ;
    t = (startTime:endTime)/fps ; %seconds
    
    [~,i1R] = min(abs(fwdFlipTimesR)) ;
    t1R = fwdFlipTimesR(i1R) ; 
    %flipTimeDiff = min(abs(fwdFlipTimesR - t1)) ;
    %[~,i1] = min(flipTimeDiff) ;
    
    if t1R <= 0 
        i2R = i1R+1 ;
        t2R = fwdFlipTimesR(i2R) ;
        phaseR = -2*pi*t1R/(t2R - t1R) ;
    elseif t1R > 0
        i2R = i1R - 1 ;
        t2R = fwdFlipTimesR(i2R) ;
        phaseR = -2*pi*t2R/(t1R - t2R) ;
    end
    
    [~,i1L] = min(abs(fwdFlipTimesL)) ;
    t1L = fwdFlipTimesL(i1L) ; 
    if t1L <= 0 
        i2L = i1L+1 ;
        t2L = fwdFlipTimesL(i2L) ;
        phaseL = -2*pi*t1L/(t2L - t1L) ;
    elseif t1L > 0
        i2L = i1L - 1 ;
        t2L = fwdFlipTimesL(i2L) ;
        phaseL = -2*pi*t2L/(t1L - t2L) ;
    end
    
    pulsePhaseArrayR(i) =phaseR ;
    pulsePhaseArrayL(i) =phaseL ;
        
end

PhaseData = table(Expr, Mov, pulsePhaseArrayR,pulsePhaseArrayL) ; 
