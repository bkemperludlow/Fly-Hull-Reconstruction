%Needs lots of commenting and cleaning up

filename = 'pitchDown_wingpitch.mat' ;
plotFlag = 1;
PitchType = 'up' ;

%For pitch up
ExprNum =  7 ; %[ 7 7 7 7 7 7 7 7 7 7 14] ;
MovNum =  {'009'} ; %{'007' '008' '009' '015' '019' '023' '024' '046' ...
%    '054' '062' '082'}; 
windowSize = 30 ; %[30 30 30 30 30 30 30 44 30 40 40]';
%Note: windowSize (defined above) fixes the window size to the left. The 
%size on the right is forced to be more or less symmetric

%For pitch down
%ExprNum = [ 7 7 7 7] ;
%MovNum =  {'006' '040' '065' '067'} ;  %Expr7Mov61 finds the wrong flip. Need to fix
%windowSize = [30 30 40 30]';
                                  
    
if (length(ExprNum) ~= length(MovNum))
    disp('Fix either ExprNum or MovNum arrays; they do not match')
    return ;
end

runSize = length(ExprNum) ;

defineConstantsScript
patchColor = [1 1 1 ] * 0.8 ; 
faceAlpha = 1 ;

pitchAcceleration = zeros(runSize,1) ;
AccelFitErr = zeros(runSize,1) ;
windowErr = zeros(runSize,1) ; 
Mov = zeros(runSize,1) ;
Expr = ExprNum' ;
%Capitalization means array. Sketchy, I know
EtaFwd = zeros(runSize,1) ;
RawEtaFwd = zeros(runSize,1) ;
EtaBack = zeros(runSize,1) ;
RawEtaBack = zeros(runSize,1) ;

for i = 1:runSize
    
    %Load data
    if strcmp(PitchType,'up')
        datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
            num2str(ExprNum(i)), '_mov_',MovNum{i}) ;
    elseif strcmp(PitchType,'down')
        datapath = strcat('F:\luca\Analysis\pitch down\Expr_',...
            num2str(ExprNum(i)), '_mov_',MovNum{i}) ;
    else
        disp('check PitchType')
        return ; 
    end
    cd(datapath)
    
    datafilename = strcat(datapath,'\Expr',num2str(ExprNum(i)), ...
        'mov',MovNum{i}, '_Data_manually_corrected.mat') ;
    
    load(datafilename) ;
    Mov(i) = str2num(MovNum{i}) ;
    
    %Check ignoreFrames and correctionTime
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
    
    if isfield(data.params,'pulseLengthMS')
        pulseLength = data.params.pulseLengthMS ;
    elseif (~isfield(data.params,'pulseLengthMS')) && (ExprNum(i) == 7)
        pulseLength = 5.8 ; %ms
    else
        pulseLength = 8 ; %ms
    end
    
    %Get flip times
    if isfield(data,'fwdFlipTimesR') && isfield(data,'backFlipTimesR') 
        fwdFlipTimesR = data.fwdFlipTimesR ; 
        backFlipTimesR = data.backFlipTimesR ; 
    else
        [fwdFlipTimesR, backFlipTimesR, fwdFlipTimesL, backFlipTimesL, ~,~, data]...
            = saveWingFlipsAndAngles(ExprNum(i),Mov(i),PitchType) ;
    end
    
    %Get the body pitch angle
    if (isfield(data, 'anglesLabFrame'))
        bodyPitch = data.anglesLabFrame(:,BETA) ;
    else
        [~, ~, ~, ~, ~, anglesLabFrame, data] = ...
            saveWingFlipsAndAngles(ExprNum(i),Mov(i),PitchType) ;
    end
    
    %exclude bad frames
    startTime = data.params.startTrackingTime ;
    endTime = data.params.endTrackingTime ;
    fps = data.params.fps ;
    t = (startTime:endTime)/fps ; %seconds
    tms = 1000*t ; %ms
    dt = 1/fps ;
    zeropoint = find(t == 0) ;
    
    bodyPitch(ignoreFrames) = NaN ; 

    ind = find(~isnan(bodyPitch)) ;
    currtvec = t(ind) ;
    currBodyPitch = bodyPitch(ind)' ;
    
    %Fit to parabola. Note maxCurrPitch refers to max pitch deviation
    earlyFrames = find(t <= .02) ; 
    if strcmp(PitchType,'up')
        [maxCurrPitch, maxCurrInd] = max(currBodyPitch(earlyFrames)) ;
    elseif strcmp(PitchType,'down')
        [maxCurrPitch, maxCurrInd] = min(currBodyPitch(earlyFrames)) ;
    end
    leftInd = maxCurrInd - windowSize(i) ;
    leftPitch = currBodyPitch(maxCurrInd - windowSize(i)) ; 
    windowTol = .5 ;
    rightInd = find(abs(currBodyPitch(maxCurrInd:end) - leftPitch) < windowTol, 1, 'first') + maxCurrInd;
    while isempty(rightInd)
        windowTol = windowTol + .1 ;
        rightInd = find(abs(currBodyPitch(maxCurrInd:end) - leftPitch) < windowTol, 1, 'first') + maxCurrInd;
    end
    polyInd = leftInd:rightInd ;
    
    bodyPitch_quadratic = fit(currtvec(polyInd)',currBodyPitch(polyInd)','poly2') ;
    bodyPitch_coeffvals = coeffvalues(bodyPitch_quadratic) ;
    bodyPitch_confInt = confint(bodyPitch_quadratic) ;
    pitchAcceleration(i) = 2*bodyPitch_coeffvals(1) ; 
    
    %Find errors due to both fitting and windowing
    AccelFitErr(i) = 2*abs(bodyPitch_confInt(1,1) - bodyPitch_coeffvals(1)) ; 
    bodyPitch_cubic = fit(currtvec(polyInd)',currBodyPitch(polyInd)','poly3') ;
    cubicCoeffvals = coeffvalues(bodyPitch_cubic) ; 
    windowErr(i) = cubicCoeffvals(1)*(windowSize(i)/8000)^3 ; 
    
    %Find changes in stroke kinematics about this parabola
    [etaFwd, etaBack, rawEtaFwd, rawEtaBack] = wingpitchAngleChange(data,bodyPitch_coeffvals,t,plotFlag) ;
    
    EtaFwd(i) = etaFwd ;
    RawEtaFwd(i) = rawEtaFwd ;
    EtaBack(i) = etaBack ;
    RawEtaBack(i) = rawEtaBack ;

    %Plot stuff
    
    tsfvec = [0 pulseLength pulseLength 0 0] ;
    
    if plotFlag == 1
        figure; 
        hold on;
        %htemp1 = plot(currtvec*1000, bodyPitch_sg, 'r-','LineWidth',2) ;
        htemp2 = plot(tms,bodyPitch,'kx','MarkerSize',4.5) ;
        ylim = get(gca,'ylim') ;
        %delete(htemp1); 
        delete(htemp2) ;
        
        avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
        hf = fill(tsfvec , avec,'y') ;
        set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
        
        %plot(currtvec*1000, bodyPitch_sg, 'r-','LineWidth',2) ;
        plot(tms,bodyPitch,'kx','MarkerSize',4.5) ;
        plot(1000*t,bodyPitch_quadratic(t),'b--','LineWidth',2)
        plot(1000*currtvec(polyInd(1)),currBodyPitch(polyInd(1)),'c.','MarkerSize',20)
        plot(1000*currtvec(polyInd(end)),currBodyPitch(polyInd(end)),'c.','MarkerSize',20)
        box on ; grid on;
        ylabel('Body Pitch Angle [deg]')
        xlabel('t [ms]')
        legend({'Pert.', '\theta_{B,raw}', 'PolyFit',...
            'Window Bound'}, 'Location', 'bestoutside') 
        title(strcat('Expr ',num2str(ExprNum(i)),' Mov ',MovNum(i)))
      
        plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
        
        set(gca,'ylim',ylim) ;
        set(gca,'xlim',xlim) ;
    end
    
end
%{
cd('F:\luca\Analysis\metadata') ;
currAcceldata = table(Expr, Mov, pitchAcceleration, AccelFitErr, windowSize, windowErr,...
    EtaFwd, EtaBack, RawEtaFwd, RawEtaBack) ;
currAcceldata.Properties.VariableUnits = {'' '' 'deg/s^2' 'deg/s^2' 'frames' 'deg/s^2'...
    'deg' 'deg' 'deg' 'deg' } ;
currAcceldata.Properties.VariableDescriptions{'Expr'} = 'Experiment Number' ;
currAcceldata.Properties.VariableDescriptions{'Mov'} = 'Movie Number' ;
currAcceldata.Properties.VariableDescriptions{'pitchAcceleration'} = '\"{\theta}' ;

if (exist(filename) == 2)
    oldAcceldata = importdata(filename) ;
    Acceldata = union(oldAcceldata,currAcceldata) ;
    save(filename,'Acceldata') ;
else
    Acceldata = currAcceldata ; 
    save(filename,'Acceldata') ;
end
%}



