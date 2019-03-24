%function [] = pitch_and_phiFront(exprNum, movNum, pertType, saveFlag, plotFlag)

defineConstantsScript

%pertType is a string. either 'Pitch Up', 'Pitch Down', 'No Perturbation',
%or 'Other'

exprNum = 13 ;
movNum  = 2 ;
saveFlag = false ;
plotFlag = true ;

%Note: pulse timing is for Janelia flies only. I think Expr6, Expr7, 
%and Expr8 all used 8 ms (pretty positive), and I know for a fact
%that later experiments used 5.8 ms
if exprNum >= 9 
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
%{
datapath = ['D:\Janelia Flies\Analysis\' pertType '\Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '\' ] ;

savefilename = [ datapath ...
    'Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '_results_temp.mat' ] ; 
%datapath = 'F:\luca\Analysis\pitch down\Expr_7_mov_067\' ;
%datafilename = [datapath 'Expr7mov67_results.mat'] ;

try 
    datafilename = [ datapath ...
    'Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '_results.mat' ] ; 
    load(datafilename) ;
catch
    datafilename = [ datapath ...
    'Expr' num2str(exprNum) 'mov' zstr num2str(movNum) '_results.mat' ] ; 
    load(datafilename) ;
end
%{
if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end
if (isfield(data,'correctionTime'))
    correctionTime=data.correctionTime;
elseif(isfield(data,'manualCorrRangeMS'))
    correctionTime=data.manualCorrRangeMS;
else 
    correctionTime = [-60 30] ;
end
%}
%}
rollEstErr = 1.5; %.5
[anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody,sp_rho, smoothed_rho, rho_t, rho_samp] ...
    = calcAngles_quick_and_dirty_mk3(data, rollEstErr, true) ;


faceAlpha  = 1 ;
a1 = -200 ; a2 = 200;
patchColor = [1 1 1 ] * 0.8;
% these vectors are used in plotting the perturbation yellow rect

t1=pulseStartMS ; t2 = pulseStartMS + pulseLengthMS ;
tsfvec = [t1 t2 t2 t1 t1 ] ;
avec   = [ a1 a1 a2 a2 a1] ;
clear t1 t2 ;

tms = t * 1000 ;
dt = 1 / data.params.fps ;
%currN = length(t);


%%
%
% [ fwdFlipIndR, backFlipIndR, fwdFlipTimesR,  backFlipTimesR,  ...
%     fwdFlipIndL,  backFlipIndL,  fwdFlipTimesL, backFlipTimesL ] = ...
%     findWingFlipTimes_mk2 (sp_phiR, sp_phiL, t) ;

phiR = -anglesBodyFrame(:, PHIR) ;

for i = 1:length(phiR)
    while phiR(i) < -20
        phiR(i) = phiR(i) + 360 ; 
    end
    while phiR(i) > 360
        phiR(i) = phiR(i) - 360 ; 
    end
end

%if (~isempty(ignoreIndR))
%    phiR(ignoreIndR) = NaN ;
%end
%}
phiL = +anglesBodyFrame(:, PHIL) ;

for i = 1:length(phiL)
    while phiL(i) < -20
        phiL(i) = phiL(i) + 360 ; 
    end
    while phiL(i) > 360
        phiL(i) = phiL(i) - 360 ; 
    end
end

%if (~isempty(ignoreIndL))
%    phiL(ignoreIndL) = NaN ;
%end

hr = figure ; plot(t, -anglesBodyFrame(:, PHIR),'ro') ;
xlabel('t') ; ylabel('\phi_R [deg]') ;
% spline smooth phiL and phiR
%ignoreIndR = unique([find(isnan(anglesBodyFrame(:, PHIR))==1)' ignoreFrames]) ;% [127 311 312 425 908 909] ; % expr7mov9

hl = figure ; plot(t, anglesBodyFrame(:, PHIL),'bo') ;
xlabel('t') ; ylabel('\phi_L [deg]') ;
% spline smooth phiL and phiR
%ignoreIndL = unique([find(isnan(anglesBodyFrame(:, PHIL))==1)'  ignoreFrames])  ; % [425 536 537 869 908 1061] ; % expr7mov9

[fwdFlipTimesR, backFlipTimesR, fwdFlipIndR, backFlipIndR, fwdFlipPhiR, backFlipPhiR, badIndicesR] = findWingFlipTimes_mk3 (t, phiR, true);
title('\phi Right') ;
[fwdFlipTimesL, backFlipTimesL, fwdFlipIndL, backFlipIndL, fwdFlipPhiL, backFlipPhiL, badIndicesL] = findWingFlipTimes_mk3 (t, phiL, true);
title('\phi Left') ;

%phiR(badIndicesR) = NaN ;
%phiL(badIndicesL) = NaN ;

anglesBodyFrame(:,PHIR) = -phiR ;
anglesBodyFrame(:,PHIL) = phiL ;

%% Smooth stroke angles

% error for pitch (used later)
ESTERR_PITCH = .6 ; %.33

% error for phi
ESTERR = 2; % measurement is (+/-)esterr about the "real" data in degrees
ESTERR_LOW = 1 ;

% smooth phiR
% use a larger tolerance (ESTERR) just for determining the flipping points.
% use a smaller tolerance for getting phi value at the flipping point

ind = ~isnan(phiR) ;
currtvec = t(ind) ;
phiR2    = phiR(ind) ;
currN = length(currtvec) ;
tol = currN * ESTERR^2 * dt ;
tol_low = currN * ESTERR_LOW^2 * dt ;
[sp_phiR, phiR_smooth] =  spaps(currtvec, phiR2, tol) ;
[sp_phiR_low, phiR_smooth_low] =  spaps(currtvec, phiR2, tol_low) ;

figure(hr)  ;
hold on ;
plot(currtvec, phiR_smooth,'k-') ;
hold off ;

% smooth phiL
ind = ~isnan(phiL) ;
currtvec = t(ind) ;
phiL2    = phiL(ind) ;
currN = length(currtvec) ;
tol = currN * ESTERR^2 * dt ;
tol_low = currN * ESTERR_LOW^2 * dt ;

[sp_phiL, phiL_smooth] =  spaps(currtvec, phiL2, tol) ;
[sp_phiL_low, phiL_smooth_low] =  spaps(currtvec, phiL2, tol_low) ;

figure(hl)  ;
hold on ;
plot(currtvec, phiL_smooth,'k-') ;
hold off ;

phiR_smooth = fnval(sp_phiR_low, t) ;
phiL_smooth = fnval(sp_phiL_low, t) ;

phiRfwdFlip = fwdFlipPhiR' ;  % fnval(sp_phiR_low, fwdFlipTimesR) ;
phiRbckFlip = backFlipPhiR' ; % fnval(sp_phiR_low, backFlipTimesR) ;
phiLfwdFlip = fwdFlipPhiL' ;  % fnval(sp_phiL_low, fwdFlipTimesL) ;
phiLbckFlip = backFlipPhiL' ; % fnval(sp_phiL_low, backFlipTimesL) ;
%%-----
close(hr) ;
close(hl) ;

%% smooth body pitch

tol = ESTERR_PITCH ;

[sp_pitch, pitch_smooth] =  mySplineSmooth(t, anglesLabFrame(:,BETA), tol) ;

%% Plot results

%plot body pitch, body roll, and front stroke angle for each movie
if plotFlag 
    hres = figure('position',[140 100 1100 800])  ;
    set(hres,'name',['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)],'numbertitle','off')
    title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
    
    s1 = subplot(2,2,1) ;
    hold on
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    plot(tms, anglesLabFrame(:,BETA),'.','color',[0.2 0.8 0.2]*0) ;
    plot(tms, fnval(sp_pitch,t)) ;
    xlabel('Time [ms]')
    ylabel('Body Pitch Angle [deg]')
    grid on ; box on ;
    set(gca,'ylim',[min(anglesLabFrame(:,BETA))-5, max(anglesLabFrame(:,BETA))+5]) ;
    set(gca,'xlim',[tms(1) tms(end)]) ;

    
    s2 = subplot(2,2,2) ;
    hold on
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    plot(rho_t*1000, rho_samp,'ko','markerfacecolor','g') ;
    plot(tms, smoothed_rho,'k-') ;
    %title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
    xlabel('Time [ms]') ;
    ylabel('Body roll angle [deg]') ;
    grid on ; box on ;
    set(gca,'ylim',[min(smoothed_rho)-5, max(smoothed_rho)+5]) ;
    set(gca,'xlim',[tms(1) tms(end)]) ;

    
    s3 = subplot(2,2,[3,4]) ;
    hold on
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    plot(fwdFlipTimesR*1000, fwdFlipPhiR, 'ro') %, 'MarkerFaceColor','r')
    plot(fwdFlipTimesL*1000, fwdFlipPhiL, 'bo') %, 'MarkerFaceColor','b')
    plot(fwdFlipTimesR*1000, mean([fwdFlipPhiR, fwdFlipPhiL],2), 'ko-', 'MarkerFaceColor','k')
    %plot(tms, phiR_smooth,'r--')
    %plot(tms, phiL_smooth,'b--')
    xlabel('Time [ms]')
    ylabel('Front Stroke Angle [deg]')
    grid on ; box on ;
    set(gca,'ylim',[min([fwdFlipPhiR; fwdFlipPhiL])-3, max([fwdFlipPhiL; fwdFlipPhiL])+3]) ;
    set(gca,'xlim',[tms(1) tms(end)]) ;

end
%% Save results

if saveFlag
    cd(datapath)
    data.fwdFlipTimesR = fwdFlipTimesR ;
    data.backFlipTimesR = backFlipTimesR ;
    data.fwdFlipTimesL = fwdFlipTimesL ;
    data.backFlipTimesL = backFlipTimesL ;
    
    data.fwdFlipPhiR = fwdFlipPhiR ;
    data.fwdFlipPhiL = fwdFlipPhiL ;
    
	data.params.pulseLengthMS = pulseLengthMS ; 
	data.anglesBodyFrame = anglesBodyFrame ;
    data.anglesLabFrame = anglesLabFrame ;
    
    data.tms = tms ;
    
    %print figure!
    savefig(hres,'Pitch and Phi_front') 
    save(savefilename,'data')
end
%
%end
%{
%% find stroke amplitude
phir_comb = [phiRfwdFlip, phiRbckFlip] ;
phir_t    = [fwdFlipTimesR, backFlipTimesR] ;
rmat = [phir_t', phir_comb'] ;
rmat = sortrows(rmat,1) ;

phir_amp = abs(diff(rmat(:,2))) ;
phir_amp_t = rmat(1:end-1,1) + diff(rmat(:,1))/2 ;
mid_stroke_r = rmat(1:end-1,2) + diff(rmat(:,2))/2 ;

phil_comb = [phiLfwdFlip, phiLbckFlip] ;
phil_t    = [fwdFlipTimesL, backFlipTimesL] ;
lmat = [phil_t', phil_comb'] ;
lmat = sortrows(lmat,1) ;

phil_amp = abs(diff(lmat(:,2))) ;
phil_amp_t = lmat(1:end-1,1) + diff(lmat(:,1))/2 ;
mid_stroke_l = lmat(1:end-1,2) + diff(lmat(:,2))/2 ;

data.phiL_amp = phil_amp ;
data.phiR_amp = phir_amp ;
save(datafilename,'data')

%%




hrl = figure('position',[140 550 1100 420])  ; hold on ;
%hf = fill(tsfvec , avec,'y') ;
%set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;

plot(tms, -anglesBodyFrame(:, PHIR),'ro') ;
plot(tms, +anglesBodyFrame(:, PHIL),'bo') ;
plot(phir_amp_t*1000, mid_stroke_r,'rs-','linewidth',2) ;
plot(phil_amp_t*1000, mid_stroke_l,'bs-','linewidth',2) ;
plot(tms, anglesLabFrame(:,BETA),'.','color',[0.2 0.8 0.2]*0) ;
plot(tms, phiR_smooth,'r-') ;
plot(tms, phiL_smooth,'b-') ;
legend({'Pert.', '\phi_R(t)','\phi_L(t)','\phi^{MID}_R','\phi^{MID}_L', '\theta_b(t)'},'location','northwest') ;
xlabel('Time [ms]') ;
ylabel('Stroke angle [deg]') ;
grid on ; box on ;
axis([tms(1) tms(end) 0 200]) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);


plot(fwdFlipTimesR*1000, phiRfwdFlip,'k^','markerfacecolor','r') ;
plot(backFlipTimesR*1000, phiRbckFlip,'kv','markerfacecolor','r') ;

plot(fwdFlipTimesL*1000, phiLfwdFlip,'k^','markerfacecolor','b') ;
plot(backFlipTimesL*1000, phiLbckFlip,'kv','markerfacecolor','b') ;

%%
hamp = figure('position',[140 100 1100 800])  ;
subplot(3,1,1) ;
hold on ;
htemp1 = plot(phir_amp_t*1000, phir_amp,'rs-','linewidth',2) ;
htemp2 = plot(phil_amp_t*1000, phil_amp,'bs-','linewidth',2) ;
ylim = get(gca,'ylim') ;
delete(htemp1) ; delete(htemp2) ;
%hf = fill(tsfvec , avec,'y') ;
%set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;

plot(phir_amp_t*1000, phir_amp,'rs-','linewidth',2) ;
plot(phil_amp_t*1000, phil_amp,'bs-','linewidth',2) ;
legend({'Pert.','\Phi_R','\Phi_L'},'location','northwest');
grid on ; box on ; ylabel('Wing stroke amplitude [deg]') ;
set(gca,'xlim',[tms(1) tms(end)]) ;
ylim(1)=max(ylim(1),100);
ylim(2)=min(ylim(2),185);
set(gca,'ylim',ylim ) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);

subplot(3,1,2) ;
htemp = plot(tms, anglesLabFrame(:,BETA),'.','color',[0.2 0.8 0.2]*0) ;
hold on ;
ylim = get(gca,'ylim') ;
set(gca,'ylim',ylim) ;
%hf = fill(tsfvec , avec,'y') ;
%set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
delete(htemp) ;
plot(tms, anglesLabFrame(:,BETA),'.','color',[0.2 0.8 0.2]*0) ;
plot(tms, pitch_smooth,'r-') ;

ylabel('Body pitch angle [deg]') ;
grid on ; box on ;
xlabel('Time [ms]') ;
set(gca,'xlim',[tms(1) tms(end)]) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);

subplot(3,1,3) ;
pitchvel = fnval( fnder(sp_pitch,1), t) ;
htemp = plot(tms, pitchvel,'-','color',[0.2 0.2 0.8]) ;
hold on ;
ylim = get(gca,'ylim') ;
set(gca,'ylim',ylim) ;
%hf = fill(tsfvec , avec*50,'y') ;
%set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
delete(htemp) ;
plot(tms, pitchvel,'-','color',[0.2 0.2 0.8],'linewidth',2) ;

ylabel('Body pitch velocity [deg/s]') ;
grid on ; box on ;
xlabel('Time [ms]') ;
set(gca,'xlim',[tms(1) tms(end)]) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);



%%
hpitch = figure ;
plot(tms, anglesLabFrame(:,BETA),'o','color',[0.2 0.8 0.2]) ;
xlabel('Time [ms]') ;
ylabel('Body pitch angle [deg]') ;
grid on ; box on ;
set(gca,'xlim',[tms(1) tms(end)]) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);

hroll = figure ;
hold on ;
plot(rho_t*1000, rho_samp,'ko','markerfacecolor','g') ;
plot(tms, smoothed_rho,'k-') ;
set(gca,'fontsize',14) ;
title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
xlabel('Time [ms]') ;
ylabel('Body roll angle [deg]') ;
grid on ; box on ;
set(gca,'xlim',[tms(1) tms(end)]) ; 
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
print(gcf,[datapath 'roll_angle'],'-dpng','-r300') ;


%axis([tms(1) tms(end) 30 70]) ;

%%
figure(hrl) ;
set([hrl hamp hpitch],'paperpositionmode','auto') ;
set(gca,'xlim',[tms(1) tms(end)]) ;
set(gca,'fontsize',14) ;
title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
print(gcf,[datapath 'phi_lr'],'-dpng','-r300') ;
set(gca,'xlim',[correctionTime(1) correctionTime(2)]) ;
print(gcf,[datapath 'phi_lr_zoom'],'-dpng','-r300') ;
set(gca,'xlim',[tms(1) tms(end)]) ;
saveas(gcf,[datapath 'phi_lr'],'fig') ;


figure(hamp) ;sb1= subplot(3,1,1) ; sb2=subplot(3,1,2) ; sb3 = subplot(3,1,3) ;
sbvec = [sb1 sb2 sb3];
set(sbvec,'xlim',[tms(1) tms(end)]) ;

set(sbvec,'fontsize',14) ;
axes(sb1) ; title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;

print(gcf,[datapath 'phi_amplitudes'],'-dpng','-r300') ;
set(sbvec,'xlim',[correctionTime(1) correctionTime(2)]) ;
print(gcf,[datapath 'phi_amplitudes_zoom'],'-dpng','-r300') ;
set(sbvec,'xlim',[tms(1) tms(end)]) ;
saveas(gcf,[datapath 'phi_amplitudes_zoom'],'fig') ;

%%
%Sam wants some graphs that will help visualize what's going on with the data.
%How about one with difference in stroke amplitude and roll angle? Then
%maybe one with stroke amplitude, body pitch, and mid stroke

if SWFlag
    
    phi_time_min = min(size(phir_amp_t,1),size(phil_amp_t,1)) ;
    phi_time_avg = mean([phir_amp_t(1:phi_time_min),phil_amp_t(1:phi_time_min)],2) ;
    
    hrollcomp = figure('position',[140 50 1100 500]) ;
    subplot(2,1,1);
    hold on ;
    htemp = plot(phi_time_avg*1000, phil_amp(1:phi_time_min)-phir_amp(1:phi_time_min),'rs-','linewidth',2) ;
    ylim = get(gca,'ylim') ;
    delete(htemp) ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    
    plot(phi_time_avg*1000, phil_amp(1:phi_time_min)-phir_amp(1:phi_time_min), 's-','Color',[0 .5 0],'linewidth',2) ;
    legend({'Pert.','\Phi_{Diff}'},'location','northwest');
    grid on ; box on ;
    set(gca,'fontsize',14) ;
    ylabel('Stroke amp. diff. [deg]') ;
    set(gca,'xlim',[correctionTime(1) correctionTime(2)]) ;
    set(gca,'ylim',ylim ) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
    
    subplot(2,1,2);
    hold on;
    plot(tms, smoothed_rho,'k-','linewidth',2) ;
    %hf = fill(tsfvec , avec,'y') ;
    %set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    set(gca,'fontsize',14) ;
    xlabel('Time [ms]') ;
    ylabel('Body roll angle [deg]') ;
    grid on ; box on ;
    set(gca,'xlim',[correctionTime(1) correctionTime(end)]) ;
    %plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    set(hrollcomp,'paperpositionmode','auto') ;
    print(gcf,[datapath 'rollcompSW'],'-dpng','-r300') ;
    saveas(gcf,[datapath 'rollcompSW'],'fig') ;
    
    hwingresponse = figure('position',[140 50 1100 900]) ;
    
    subplot(3,1,1);
    hold on;
    set(gca,'fontsize',14) ;
    htemp1 = plot(phir_amp_t*1000, phir_amp,'rs-','linewidth',2) ;
    htemp2 = plot(phil_amp_t*1000, phil_amp,'bs-','linewidth',2) ;
    ylim = get(gca,'ylim') ;
    delete(htemp1) ; delete(htemp2) ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    
    plot(phir_amp_t*1000, phir_amp,'rs-','linewidth',2) ;
    plot(phil_amp_t*1000, phil_amp,'bs-','linewidth',2) ;
    legend({'Pert.','\Phi_R','\Phi_L'},'location','northwest');
    grid on ; box on ; ylabel('Wing stroke amplitude [deg]') ;
    set(gca,'xlim',[correctionTime(1) correctionTime(2)]) ;
    ylim(1)=max(ylim(1),100);
    ylim(2)=min(ylim(2),185);
    set(gca,'ylim',ylim ) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
    
    subplot(3,1,2);
    hold on ;
    set(gca,'fontsize',14) ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    plot(phir_amp_t*1000, mid_stroke_r,'rs--','linewidth',1) ;
    plot(phil_amp_t*1000, mid_stroke_l,'bs--','linewidth',1) ;
    plot(phi_time_avg*1000, mean([mid_stroke_l(1:phi_time_min), mid_stroke_r(1:phi_time_min)],2),'ks-','linewidth',2) ;
    legend({'Pert.','\Phi_R Mid','\Phi_L Mid','Avg. Mid'},'location','northwest');
    xlabel('Time [ms]') ;
    ylabel('Midstroke angle [deg]') ;
    grid on ; box on ;
    axis([correctionTime(1) correctionTime(2) 80 115]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    
    
    subplot(3,1,3);
    set(gca,'fontsize',14) ;
    htemp = plot(tms, anglesLabFrame(:,BETA),'.','color',[0.2 0.8 0.2]*0) ;
    hold on ;
    ylim = get(gca,'ylim') ;
    set(gca,'ylim',ylim) ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    delete(htemp) ;
    plot(tms, anglesLabFrame(:,BETA),'.','color',[0.2 0.8 0.2]*0) ;
    plot(tms, pitch_smooth,'r-') ;
    
    ylabel('Body pitch angle [deg]') ;
    grid on ; box on ;
    xlabel('Time [ms]') ;
    set(gca,'xlim',[correctionTime(1) correctionTime(2)]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    
    set(hwingresponse,'paperpositionmode','auto') ;
    print(gcf,[datapath 'wing_responseSW'],'-dpng','-r300') ;
    saveas(gcf,[datapath 'wing_responseSW'],'fig') ;
end
%}