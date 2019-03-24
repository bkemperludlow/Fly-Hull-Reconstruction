function controller_fit_struct = ...
    controllerFit_v3(data, ExprNum, MovNum, flyType, x_0, lb, ub, debugFlag, plotFlag) 

% x = input vector of the things I'm fitting for. Should
%       be of the form [K_i, K_p, deltaT, K] 
%{
ExprNum = 3 ;
MovNum = 37 ;

if MovNum < 10
    zstr = '00' ;
elseif MovNum < 100 ;
    zstr = '0' ;
else
    zstr = '' ;
end

movieStr = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum)] ;
dataPath = ['G:\Janelia Flies\kir2.1 flies\Analysis\Pitch Up\' movieStr] ;
%dataPath = ['F:\luca\Analysis\pitch up\' movieStr] ;
cd(dataPath)
fileName = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected.mat'] ;
load(fileName) ;
%}

if nargin < 2
    ExprNum = nan ;
    MovNum = nan ;
    flyType = 'unspecified' ;
    plotFlag = true ;
    debugFlag = true ;
end

controller_fit_struct = struct('ExprNum',[],'MovNum',[],'K_i',[],'K_p',[],...
    'deltaT',[],'K',[],'rms',[],'c_pitch',[],'sp_phiR',[],...
    'sp_phiL',[],'phiEstErr',[],'deltaPhiFront',[],'fwdFlipTimes',[],'backFlipTimes', [],...
    't',[],'flyType',[]) ;

if flyType == 1
    plotColor = [.7 0 0 ] ;
    flyTypeStr = 'experimental' ;
elseif flyType == 2
    plotColor = [0 .7 0 ] ;
    flyTypeStr = 'control' ;
else
    plotColor = .5*[1 1 1] ;
end

defineConstantsScript
patchColor = [1 1 1 ] * 0.8;

if isfield(data,'manualCorrRangeMS')
    %manualCorrRangeMS = data.manualCorrRangeMS ;
    manualCorrRangeMS_start = max([data.manualCorrRangeMS(1) -10]) ;
    manualCorrRangeMS_end = min([data.manualCorrRangeMS(2) 40]) ;
    manualCorrRangeMS = [manualCorrRangeMS_start, manualCorrRangeMS_end ] ;
    %manualCorrRangeMS = [-20 55] ;
    %manualCorrRangeMS = [-10 55] ;
else
    manualCorrRangeMS = [-10 30] ;
end
manualCorrRange = manualCorrRangeMS / 1000 ; 

fwdFlipTimesR = data.fwdFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
backFlipTimesR = data.backFlipTimesR ;
backFlipTimesL = data.backFlipTimesL ;


correctedIndR = find(fwdFlipTimesR > manualCorrRange(1) & fwdFlipTimesR < manualCorrRange(2)) ;
correctedIndL = find(fwdFlipTimesL > manualCorrRange(1) & fwdFlipTimesL < manualCorrRange(2)) ;
correctedIndBackR = find(backFlipTimesR > manualCorrRange(1) & backFlipTimesR < manualCorrRange(2)) ;
correctedIndBackL = find(backFlipTimesL > manualCorrRange(1) & backFlipTimesL < manualCorrRange(2)) ;

fwdFlipTimesR = fwdFlipTimesR(correctedIndR) ;
fwdFlipTimesL = fwdFlipTimesL(correctedIndL) ;
backFlipTimesR = backFlipTimesR(correctedIndBackR) ;
backFlipTimesL = backFlipTimesL(correctedIndBackL) ;
%fwdFlipTimesL = fwdFlipTimesR ; 
%fwdFlipTimesR = fwdFlipTimesL ; 

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

%t = t(1:end-10) ;
%phiR = phiR(1:end-10)  ;
%phiL = phiL(1:end-10)   ;
%bodyPitch = bodyPitch(1:end-10) ;

phiEstErr = 1 ;
[sp_phiR, ~, ~] = mySplineSmooth(t(~isnan(phiR)),phiR(~isnan(phiR)),phiEstErr) ;
[sp_phiL, ~, ~] = mySplineSmooth(t(~isnan(phiL)),phiL(~isnan(phiL)),phiEstErr) ;


%low-pass butterworth filter for data
% d1 = designfilt('lowpassiir','FilterOrder',8,'SampleRate',8000, ...
%     'HalfPowerFrequency',150,'DesignMethod','butter'); %hpf = 100
d1 = designfilt('lowpassiir','FilterOrder',3,'SampleRate',8000, ...
    'HalfPowerFrequency',100,'DesignMethod','butter'); %hpf = 100
pitch_filt = filtfilt(d1,bodyPitch) ;
c_pitch = fit(t',pitch_filt,'cubicinterp');
pitch_smoothed = c_pitch(t) ;
pitchVel_smoothed = differentiate(c_pitch, t) ; 
c_pitchVel = fit(t',pitchVel_smoothed,'linearinterp');

%pitchEstErr = .48 ; %.5
%[sp_pitch, ~,~] = mySplineSmooth(t,bodyPitch,pitchEstErr) ;
%pitch_smoothed = fnval(sp_pitch, t) ;

if ~isfield(data,'oneWing')
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
elseif strcmp(data.oneWing,'L')
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
    phiFrontR = phiFrontL ;
    fwdFlipTimesR = fwdFlipTimesL ; 
elseif strcmp(data.oneWing,'R')
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = phiFrontR ;
    fwdFlipTimesL = fwdFlipTimesR ; 
else 
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
end
    
    

if length(fwdFlipTimesR) == length(fwdFlipTimesL)
    phiFront = (phiFrontR + phiFrontL) / 2 ;
    phiFront_errorBar = (phiFrontR - phiFrontL) / 2 ;
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL ) /2 ;
elseif length(fwdFlipTimesR) < length(fwdFlipTimesL)
    idx = zeros(length(fwdFlipTimesR),1) ;
    for q = 1:length(fwdFlipTimesR)
        [~,minInd] = min(abs(fwdFlipTimesL - fwdFlipTimesR(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL(idx)) / 2 ;
    phiFront = (phiFrontR + phiFrontL(idx)) / 2 ;
    phiFront_errorBar = (phiFrontR + phiFrontL(idx)) / 2 ;
elseif length(fwdFlipTimesL) < length(fwdFlipTimesR)
    idx = zeros(length(fwdFlipTimesL),1) ;
    for q = 1:length(fwdFlipTimesL)
        [~,minInd] = min(abs(fwdFlipTimesR - fwdFlipTimesL(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesL + fwdFlipTimesR(idx)) / 2 ;
    phiFront = (phiFrontL + phiFrontR(idx)) / 2 ;
    phiFront_errorBar = (phiFrontL - phiFrontR(idx)) / 2 ;
end

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
    %phiBack = (phiBackL + phiBackR(idx)) / 2 ;
end


if debugFlag
    figure ; hold on
    plot(fwdFlipTimesR*1000, phiFrontR,'ro')
    plot(fwdFlipTimesL*1000, phiFrontL,'bo')
    plot(fwdFlipTimes*1000, phiFront ,'ko', 'markerfacecolor', 'k')
    plot(t*1000, fnval(sp_phiL,t),'b-')
    plot(t*1000, fnval(sp_phiR,t),'r-')
    %plot(t*1000, phiR, 'rsq')
    %plot(t*1000, phiL, 'bsq')
    xlabel('Time [ms]')
    ylabel('Front Stroke Angle, \Phi , [deg]')
    
    axis tight
    
    figure ; hold on
    plot(t*1000, bodyPitch, 'k.')
    plot(t*1000, pitch_smoothed, 'r-')
    plot(t*1000, pitch_filt, 'b-')
    xlabel('Time [ms]')
    ylabel('Body Pitch, \theta_b , [deg]')
    axis tight
    
    %keyboard ;
end

prePertFlipInd = find( fwdFlipTimes < 0 ) ;
phiFrontAvg = mean(phiFront(prePertFlipInd)) ;
phiFront_meansub = phiFront - phiFrontAvg ; 

%as intial guesses, use more or less what's in the pitch paper
%K_i_guess = 0.0 ;%0.3 ; %0.3 ; %unitless
%K_p_guess = 0.008 ; %seconds
%deltaT_guess = 0.006 ; %seconds
%K_guess = 0 ; %degrees
%theta_0_guess = 45 ; 
%deltaT_2_guess = 0.03; %seconds

%x_0 = [K_i_guess K_p_guess deltaT_guess K_guess deltaT_2_guess]; %initial guess for x
%x_0 = [K_i_guess K_p_guess deltaT_guess K_guess ]; %initial guess for x

%lb = [0 0 0 -7 0 ] ;
%ub = [1 .01 .05 7 .05 ] ;
%lb = [-2 0 .002 -5 ] ;
%ub = [2 .1 .025 5 ] ;

%controllerEq = @(K_i, K_p, deltaT, K, phiFront_meansub,fwdFlipTimes,c_pitch, c_pitchVel) ...
%    K + K_p*(c_pitchVel(t-deltaT))+ K_i*(c_pitch(t-deltaT)- c_pitch(0)) ;
options = optimoptions(@fmincon,'Algorithm','sqp','Display','off') ; %sqp 'interior-point'
%[x, fval] = fmincon(@(x)controller_residuals(x,phiFront_meansub,fwdFlipTimes,sp_pitch),x_0,[],[],[],[],lb,ub,[],options) ;
[x, fval] = fmincon(@(x)controller_residuals_v2(x,phiFront_meansub,fwdFlipTimes,c_pitch),x_0,[],[],[],[],lb,ub,[],options) ;


K_i = x(1) ;
K_p = x(2) ;
deltaT = x(3) ; 
K = x(4) ; 
%deltaT_2 = x(5)  
%theta_0 = x(5) 

rms = sqrt(fval/length(fwdFlipTimes)) ;

%controlPred = K + K_p*fnval(fnder(sp_pitch,1),t-deltaT)...
%        + K_i*(fnval(sp_pitch,t-deltaT)- fnval(sp_pitch,0)) ;
pitchVel = differentiate(c_pitch,t-deltaT) ;
controlPred = K + K_p*pitchVel...
        + K_i*(c_pitch(t-deltaT)- c_pitch(0)) ;
    
proportionalTerm = K_p*pitchVel ; 
integralTerm = K_i*(c_pitch(t-deltaT)- c_pitch(0)) ; 


%controlPred_atFlips = K + K_p*fnval(fnder(sp_pitch,1),fwdFlipTimes-deltaT)...
%        + K_i*(fnval(sp_pitch,fwdFlipTimes-deltaT)- fnval(sp_pitch,0)) ;

t_range = manualCorrRangeMS(1) : 0.125 : manualCorrRangeMS(2) ;

controller_fit_struct.ExprNum = ExprNum ;
controller_fit_struct.MovNum = MovNum ;
controller_fit_struct.K_i = K_i ;
controller_fit_struct.K_p = K_p ;
controller_fit_struct.deltaT = deltaT ;
controller_fit_struct.K = K ;
controller_fit_struct.rms = rms ;
controller_fit_struct.flyType = flyTypeStr ;
controller_fit_struct.c_pitch = c_pitch ;
controller_fit_struct.sp_phiR = sp_phiR ;
controller_fit_struct.sp_phiL = sp_phiL ;
%controller_fit_struct.pitchEstErr = pitchEstErr ;
controller_fit_struct.phiEstErr = phiEstErr ;
controller_fit_struct.fwdFlipTimes = fwdFlipTimes ;
controller_fit_struct.backFlipTimes = backFlipTimes ;
controller_fit_struct.deltaPhiFront = phiFront_meansub ;
controller_fit_struct.t = t_range ;


if plotFlag 
    ylim = [min(phiFront_meansub)-3 , max(phiFront_meansub)+3] ;
    tsfvec = [0 7 7 0 0] ;
    hcontrib = figure ;
    hold on
    %set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'Position', [500 500 420 140]);
    set(gcf, 'Color', 'w')
    set(gcf, 'Renderer','opengl')
    set(gcf,'PaperPositionMode','auto')
    %set(hcontrib,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' Contributions'],'numbertitle','off')
    
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',1,'edgecolor','none') ;
    set(hf,'HandleVisibility','off')
    
    plot(1000*t,controlPred,'Color',plotColor,'LineWidth',2.5)
    plot(1000*t,proportionalTerm,'Color',0.6*[1 1 1],'LineWidth',1.5)
    plot(1000*t,integralTerm,'k--','LineWidth',1.5)
    %shadedErrorBar(t*1000, controlPred, 2*ones(size(controlPred)),{'-','LineWidth',2.5,'Color',plotColor},1) ;
    %errorbar(1000*fwdFlipTimes, phiFront_meansub, 2*ones(size(phiFront_meansub)), 'ko','markerfacecolor','k') ;
    errorbar(1000*fwdFlipTimes, phiFront_meansub, phiFront_errorBar,'ko','markerfacecolor','k') ;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    axis tight ;
    set(gca, 'xlim', 1000*[manualCorrRange(1), max(fwdFlipTimes)])
    set(gca,'ylim',ylim)
    xlabel('Time [ms]')
    ylabel('\Delta \Phi_{front} [deg]')
    %legend({'PI','P','I','Data'})
end

end 