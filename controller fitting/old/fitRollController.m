function [controller_fit_struct, h_contrib, h_debug] = ...
    fitRollController(data, ExprNum, MovNum, flyType, debugFlag, plotFlag) 

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
    'deltaT',[],'rms',[],'sp_rho',[],'rhoEstErr',[],'phiAmpDiff', [], 'phiAmpR',[],...
    'phiAmpL',[],'phiAmpTimesR',[],'phiAmpTimesL',[],'phiAmpTimes',[],'t',[],...
    'flyType',[]) ;

if flyType == 1
    %plotColor = [.7 0 0 ] ;
    plotColor = [191 87 0]/255 ;
    flyTypeStr = 'experimental' ;
elseif flyType == 2
    %plotColor = [0 .7 0 ] ;
    plotColor = [0 119 204]/255 ;
    flyTypeStr = 'control' ;
else
    plotColor = .5*[1 1 1] ;
end

defineConstantsScript
patchColor = [1 1 1 ] * 0.8;


if isfield(data,'manualCorrRangeMS')
    %manualCorrRangeMS = data.manualCorrRangeMS ;
    manualCorrRangeMS_start = max([data.manualCorrRangeMS(1) -10]) ;
    manualCorrRangeMS_end = min([data.manualCorrRangeMS(2) 50]) ;
    manualCorrRangeMS = [manualCorrRangeMS_start, manualCorrRangeMS_end ] ;
   %manualCorrRangeMS = [-20 55] ;
    %manualCorrRangeMS = [-10 30] ;
else
    manualCorrRangeMS = [-10 50] ;
end
manualCorrRange = manualCorrRangeMS / 1000 ; 

phir_amp_t = data.phir_amp_t ;
phil_amp_t = data.phil_amp_t ;
phiR_amp = data.phiR_amp ;
phiL_amp = data.phiL_amp ;
bodyRoll = data.anglesLabFrame(:,RHO) ;
fwdFlipTimesR = data.fwdFlipTimesR ; 
backFlipTimesR = data.backFlipTimesR ; 

correctedIndR = find(phir_amp_t > manualCorrRange(1) & phir_amp_t < manualCorrRange(2)) ;
correctedIndL = find(phil_amp_t > manualCorrRange(1) & phil_amp_t < manualCorrRange(2)) ;

phir_amp_t = phir_amp_t(correctedIndR) ;
phil_amp_t = phil_amp_t(correctedIndL) ;
phiR_amp = phiR_amp(correctedIndR) ;
phiL_amp = phiL_amp(correctedIndL) ;


t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;

rhoEstErr = 0.01 ; %it's already the result of a spline
[sp_rho, ~, ~] = mySplineSmooth(t,bodyRoll,rhoEstErr) ;


if isfield(data,'oneWing') && strcmp(data.oneWing,'L')
    phiR_amp = phiL_amp ; 
    phir_amp_t = phil_amp_t ;
elseif isfield(data,'oneWing') && strcmp(data.oneWing,'R')
    phiL_amp = phiR_amp ; 
    phil_amp_t = phir_amp_t ;
end
    
    
if length(phir_amp_t) == length(phil_amp_t)
    phiAmpDiff = phiR_amp - phiL_amp ;
    phiAmpTimes = (phil_amp_t + phir_amp_t ) /2 ;
elseif length(phir_amp_t) < length(phil_amp_t)
    idx = zeros(length(phir_amp_t),1) ;
    for q = 1:length(phir_amp_t)
        [~,minInd] = min(abs(phil_amp_t - phir_amp_t(q))) ; 
        idx(q) = minInd ;
    end
    phiAmpTimes = (phir_amp_t + phil_amp_t(idx)) / 2 ;
    phiAmpDiff = (phiR_amp - phiL_amp(idx)) / 2 ;
elseif length(phil_amp_t) < length(phir_amp_t)
    idx = zeros(length(phil_amp_t),1) ;
    for q = 1:length(phil_amp_t)
        [~,minInd] = min(abs(phir_amp_t - phil_amp_t(q))) ; 
        idx(q) = minInd ;
    end
    phiAmpTimes = (phir_amp_t(idx) + phil_amp_t) / 2 ;
    phiAmpDiff = (phiR_amp(idx) - phiL_amp) / 2 ;
end


%as intial guesses, use more or less what's in the roll paper
K_i_guess = 0.7 ;%0.3 ; %0.3 ; %unitless
K_p_guess = 0.006 ; %seconds
deltaT_guess = 0.006 ; %seconds


%x_0 = [K_i_guess K_p_guess deltaT_guess K_guess deltaT_2_guess]; %initial guess for x
x_0 = [K_i_guess K_p_guess deltaT_guess]; %initial guess for x

%lb = [0 0 0 -7 0 ] ;
%ub = [1 .01 .05 7 .05 ] ;
lb = [-2 0 .0035] ;
ub = [2 .1 .015] ;

options = optimoptions(@fmincon,'Algorithm','sqp') ; %sqp 'interior-point'
%[x, fval] = fmincon(@(x)controller_residuals(x,phiFront_meansub,fwdFlipTimes,sp_pitch),x_0,[],[],[],[],lb,ub,[],options) ;
[x, fval] = fmincon(@(x)controller_residuals_roll(x,phiAmpDiff,phiAmpTimes,sp_rho),x_0,[],[],[],[],lb,ub,[],options) ;

K_i = x(1) 
K_p = x(2) 
deltaT = x(3) 
%K = x(4) 
%deltaT_2 = x(5)  
%theta_0 = x(5) 

rms = sqrt(fval/length(phiAmpTimes)) 

%controlPred = K + K_p*fnval(fnder(sp_pitch,1),t-deltaT)...
%        + K_i*(fnval(sp_pitch,t-deltaT)- fnval(sp_pitch,0)) ;
rollVel = fnval(fnder(sp_rho,1),t-deltaT) ;
controlPred = K_p*rollVel + K_i*fnval(sp_rho, t - deltaT) ;
    
proportionalTerm = K_p*rollVel ; 
integralTerm = K_i*fnval(sp_rho, t - deltaT) ; 


%controlPred_atFlips = K + K_p*fnval(fnder(sp_pitch,1),fwdFlipTimes-deltaT)...
%        + K_i*(fnval(sp_pitch,fwdFlipTimes-deltaT)- fnval(sp_pitch,0)) ;

t_range = manualCorrRangeMS(1) : 0.125 : manualCorrRangeMS(2) ;

controller_fit_struct.ExprNum = ExprNum ;
controller_fit_struct.MovNum = MovNum ;
controller_fit_struct.K_i = K_i ;
controller_fit_struct.K_p = K_p ;
controller_fit_struct.deltaT = deltaT ;
%controller_fit_struct.K = K ;
controller_fit_struct.rms = rms ;
controller_fit_struct.flyType = flyTypeStr ;
controller_fit_struct.sp_rho = sp_rho ;
%controller_fit_struct.sp_phiR = sp_phiR ;
%controller_fit_struct.sp_phiL = sp_phiL ;
%controller_fit_struct.pitchEstErr = pitchEstErr ;
%controller_fit_struct.phiEstErr = phiEstErr ;
controller_fit_struct.phiAmpTimes = phiAmpTimes ;
controller_fit_struct.phiAmpDiff = phiAmpDiff ;
controller_fit_struct.t = t_range ;

if debugFlag
    xlim = 1000*[manualCorrRange(1), max(phiAmpTimes)] ; 
    h_debug = figure('Position',[680, 390, 560, 588],'PaperPositionMode','auto') ; 
    
    subplot(4,1,1)
    hold on
    plot(phir_amp_t*1000, phiR_amp,'ro-')
    plot(phil_amp_t*1000, phiL_amp,'bo-')
    set(gca,'xlim',xlim)
    ylabel({'Stroke Amp.', '\Phi , [deg]'})
    legend({'\Phi_R','\Phi_L'},'location', 'northwest')
    grid on
    
    subplot(4,1,2)
    plot(phiAmpTimes*1000, phiAmpDiff ,'ko-', 'markerfacecolor', 'k')
    %xlabel('Time [ms]')
    ylabel({'Stroke Amp. Diff.', '\Delta\Phi , [deg]'})
    set(gca,'xlim',xlim)
    grid on

    subplot(4,1,3)
    hold on
    plot(t*1000, bodyRoll, 'k.')
    plot(t*1000, fnval(sp_rho, t), 'r-')
    %xlabel('Time [ms]')
    ylabel({'Body Roll', '\rho , [deg]'})
    set(gca,'xlim',xlim)
    %set(gca,'xlim',[phiAmpTimes(1)*1000  phiAmpTimes(end)*1000])
    max_rho = max(bodyRoll((t > phiAmpTimes(1) & t < phiAmpTimes(end)))) ; 
    min_rho = min(bodyRoll((t > phiAmpTimes(1) & t < phiAmpTimes(end)))) ; 
    set(gca,'ylim',[(min_rho-5) (max_rho+5)]) ; 
    grid on
    %keyboard ;
    
    sp_4 = subplot(4,1,4) ;
    ylim = [min(phiAmpDiff)-10 , max(phiAmpDiff)+10] ;
    tsfvec = [0 7 7 0 0] ;
    hold on
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',.7) ;
    set(hf,'HandleVisibility','off')
    
    h_pi = plot(1000*t,controlPred,'Color',plotColor,'LineWidth',2.5) ;
    h_p = plot(1000*t,proportionalTerm,'Color',0.6*[1 1 1],'LineWidth',1.5) ;
    h_i = plot(1000*t,integralTerm,'k--','LineWidth',1.5) ;
    %shadedErrorBar(t*1000, controlPred, 2*ones(size(controlPred)),{'-','LineWidth',2.5,'Color',plotColor},1) ;
    %errorbar(1000*fwdFlipTimes, phiFront_meansub, 2*ones(size(phiFront_meansub)), 'ko','markerfacecolor','k') ;
    h_data = plot(1000*phiAmpTimes, phiAmpDiff, 'ko','markerfacecolor','k') ;
    axis tight ;
    set(gca, 'xlim', xlim)
    set(gca,'ylim',ylim)
    legend([h_pi, h_p, h_i, h_data], {'PI', 'Prop.', 'Int.', 'Data'},...
        'location','northoutside','orientation','horizontal')
    xlabel('Time [ms]')
    ylabel('\Delta \Phi [deg]')
    grid on
else 
    h_debug = nan ;
end


if plotFlag 
    ylim = [min(phiAmpDiff)-10 , max(phiAmpDiff)+10] ;
    tsfvec = [0 7 7 0 0] ;
    h_contrib = figure ;
    hold on
    %set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'Position', [500 500 420 140]);
    %set(gcf, 'Color', 'w')
    %set(gcf, 'Renderer','opengl')
    set(gcf,'PaperPositionMode','auto')
    %set(hcontrib,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' Contributions'],'numbertitle','off')
    
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',1,'edgecolor','none') ;
    set(hf,'HandleVisibility','off')
    
    plot(1000*t,controlPred,'Color',plotColor,'LineWidth',3)
    plot(1000*t,proportionalTerm,'Color',0.6*[1 1 1],'LineWidth',1.5)
    plot(1000*t,integralTerm,'k--','LineWidth',1.5)
    %shadedErrorBar(t*1000, controlPred, 2*ones(size(controlPred)),{'-','LineWidth',2.5,'Color',plotColor},1) ;
    errorbar(1000*phiAmpTimes, phiAmpDiff, 4*ones(size(phiAmpTimes)), 'ko', 'markerfacecolor','k') ;
    %plot(1000*phiAmpTimes, phiAmpDiff, 'ko','markerfacecolor','k') ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
    axis tight ;
    set(gca, 'xlim', 1000*[manualCorrRange(1), max(phiAmpTimes)])
    set(gca,'ylim',ylim)
    xlabel('Time [ms]')
    ylabel('\Delta \Phi [deg]')
    %legend({'PI','P','I','Data'})
else
    h_contrib = nan ; 
end

end 