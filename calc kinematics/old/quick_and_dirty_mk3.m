
%{
load results

defineConstantsScript;
%}

%{
clear data ;
for movNum = [7 10] ;
    quick_and_dirty ;
    close all ;
end
%}

SWFlag = 0 ;
SAVE_FLAG = true ; 
largePertFlag = true ; 

defineConstantsScript

%{
if exprNum == 7
    pulseLengthMS = 5.8 ;
else 
    pulseLengthMS = 8 ;
end
%}
%filterSpec = 'I:\B1 Yaw Data\Yaw analysis\Saumya Manually Corrected\' ; 
filterSpec = 'I:\Fly Data\VNC Motor Lines\' ; 
%filterSpec = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ; 

filepath_cell = uipickfiles('FilterSpec',filterSpec,'Type',{'*.mat', 'MAT-files' },'NumFiles',1) ; 
filepath = filepath_cell{1} ; 

filepath_split = strsplit(filepath,'\') ; 
datafilename = filepath_split{end} ;
datapath = filepath(1:end-length(datafilename)) ; 

if contains(datafilename,'Expr_')
    datafilename_split = strsplit(datafilename,'_') ;
    exprNum = str2double(datafilename_split{2});
    movNum = str2double(datafilename_split{4}) ;
else
    expr_idx = strfind(datafilename,'Expr') ; 
    mov_idx = strfind(datafilename,'mov') ; 
    
    exprNum = str2double(datafilename((expr_idx+4):(mov_idx-1))) ; 
    movNum = str2double(datafilename((mov_idx+3):(mov_idx+5))) ; 
end

if (movNum<10)
    zstr = '00' ;
elseif (movNum<100)
    zstr = '0' ;
else
    zstr = '';
end

pulseLengthMS = 7 ;
pulseStartMS = 0 ;

data = importdata(filepath) ; 
%
%datapath = ['F:\Sam\Janelia Flies\Analysis\No perturbation\Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '\' ] ;

%datafilename = [ datapath ...
%    'Expr' num2str(exprNum) 'mov' zstr num2str(movNum) '_Data_manually_corrected.mat' ] ; 

%datapath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Up\Expr_11_mov_073\';
%datafilename = [datapath 'Expr11mov073_Data_manually_corrected.mat'] ;
%datafilename = [datapath 'Expr_7_mov_020_test.mat'] ;
%datafilename = [datapath 'Expr7mov020_Data_manually_corrected_onlyPhi.mat'] ;

%cd(datapath)
%load(datafilename) ;
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
    correctionTime = [-49 53] ;
end
rollEstErr = .05 ; %.3
% [anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody,sp_rho, smoothed_rho, rho_t, rho_samp] ...
%     = calcAngles_quick_and_dirty_mk2_C7(data, rollEstErr, false, bigPertFlag) ;
[anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody, sp_rho,...
    smoothed_rho, rho_t, rho_samp] = ...
    calcAnglesRaw_Sam(data, false,largePertFlag) ; 

faceAlpha  = 1 ;
a1 = -200 ; a2 = 200;
patchColor = [1 1 1 ] * 0.8;
% these vectors are used in plotting the perturbation yellow rect
%{
t1=pulseStartMS ; t2 = pulseStartMS + pulseLengthMS ;
tsfvec = [t1 t2 t2 t1 t1 ] ;
avec   = [ a1 a1 a2 a2 a1] ;
clear t1 t2 ;
%}

%%
%
% [ fwdFlipIndR, backFlipIndR, fwdFlipTimesR,  backFlipTimesR,  ...
%     fwdFlipIndL,  backFlipIndL,  fwdFlipTimesL, backFlipTimesL ] = ...
%     findWingFlipTimes_mk2 (sp_phiR, sp_phiL, t) ;
hr = figure ; plot(t, -anglesBodyFrame(:, PHIR),'ro') ;
xlabel('t') ; ylabel('\phi_R [deg]') ;
% spline smooth phiL and phiR
ignoreIndR = unique([find(isnan(anglesBodyFrame(:, PHIR))==1)' ignoreFrames]) ;% [127 311 312 425 908 909] ; % expr7mov9

hl = figure ; plot(t, anglesBodyFrame(:, PHIL),'bo') ;
xlabel('t') ; ylabel('\phi_L [deg]') ;
% spline smooth phiL and phiR
ignoreIndL = unique([find(isnan(anglesBodyFrame(:, PHIL))==1)'  ignoreFrames])  ; % [425 536 537 869 908 1061] ; % expr7mov9

phiR = -anglesBodyFrame(:, PHIR) ;
%phiR = phiR + 360 ;

for i = 1:length(phiR)
    while phiR(i) < 0 
        phiR(i) = phiR(i) + 360 ; 
    end
    while phiR(i) > 270
        phiR(i) = phiR(i) - 360 ; 
    end
end

if (~isempty(ignoreIndR))
    phiR(ignoreIndR) = NaN ;
end

phiL = +anglesBodyFrame(:, PHIL) ;

for i = 1:length(phiL)
    while phiL(i) < -90
        phiL(i) = phiL(i) + 360 ; 
    end
    while phiL(i) > 270
        phiL(i) = phiL(i) - 360 ; 
    end
end

if (~isempty(ignoreIndL))
    phiL(ignoreIndL) = NaN ;
end

[~, hampelR] = hampel(phiR, 5) ;
[~, hampelL] = hampel(phiL, 5) ;

phiR(hampelR) = NaN ;
phiL(hampelL) = NaN ; 

[fwdFlipTimesR, backFlipTimesR, fwdFlipIndR, backFlipIndR, fwdFlipPhiR, backFlipPhiR, badIndicesR] = findWingFlipTimes_mk3 (t, phiR, true);
title('\phi Right') ;
[fwdFlipTimesL, backFlipTimesL, fwdFlipIndL, backFlipIndL, fwdFlipPhiL, backFlipPhiL, badIndicesL] = findWingFlipTimes_mk3 (t, phiL, true);
title('\phi Left') ;

phiR(badIndicesR) = NaN ;
phiL(badIndicesL) = NaN ;

anglesBodyFrame(:,PHIR) = -phiR ;
anglesBodyFrame(:,PHIL) = phiL ;

data.fwdFlipTimesR = fwdFlipTimesR ;
data.fwdFlipIndR = fwdFlipIndR ; 
data.backFlipTimesR = backFlipTimesR ;
data.backFlipIndR = backFlipIndR ; 
data.fwdFlipTimesL = fwdFlipTimesL ;
data.fwdFlipIndL = fwdFlipIndL ; 
data.backFlipTimesL = backFlipTimesL ;
data.backFlipIndL = backFlipIndL ; 
%data.params.pulseLengthMS = pulseLengthMS ; 
data.anglesBodyFrame = anglesBodyFrame ;
data.anglesLabFrame = anglesLabFrame ;

%save(datafilename,'data')
%%



dt = 1/8000 ;

% error for pitch (used later)
ESTERR_PITCH = .33 ;

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


%%-----
close(hr) ;
close(hl) ;

%% smooth body pitch

dt = 1 / data.params.fps ;
currN = length(t);

tol = currN * ESTERR_PITCH^2 * dt ;

[sp_pitch, pitch_smooth] =  spaps(t, anglesLabFrame(:,BETA), tol) ;

phiRfwdFlip = fwdFlipPhiR' ;  % fnval(sp_phiR_low, fwdFlipTimesR) ;
phiRbckFlip = backFlipPhiR' ; % fnval(sp_phiR_low, backFlipTimesR) ;
phiLfwdFlip = fwdFlipPhiL' ;  % fnval(sp_phiL_low, fwdFlipTimesL) ;
phiLbckFlip = backFlipPhiL' ; % fnval(sp_phiL_low, backFlipTimesL) ;

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

data.phil_amp_t = phil_amp_t ;
data.phir_amp_t = phir_amp_t ;
data.phiL_amp = phil_amp ;
data.phiR_amp = phir_amp ;

if SAVE_FLAG
    save(filepath,'data')
end

%%
phiR_smooth = fnval(sp_phiR_low, t) ;
phiL_smooth = fnval(sp_phiL_low, t) ;
tms = t * 1000 ;



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
%legend({'\phi_R(t)','\phi_L(t)','\phi^{MID}_R','\phi^{MID}_L', '\theta_b(t)'},'location','northwest') ;
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

[anglesMat_R, smooth_anglesMat_R, sp_phi_R, sp_theta_R, sp_psi_R ] = ... 
    smoothWingAngles(data, 'R') ; 
[anglesMat_L, smooth_anglesMat_L, sp_phi_L, sp_theta_L, sp_psi_L ] = ... 
    smoothWingAngles(data, 'L') ; 

h_wa = figure('position',[140 100 1100 800])  ;
title('Wing Angles')
%----------------
ax1 = subplot(3,1,1) ;
hold on ;
plot(tms, -anglesBodyFrame(:, PHIR),'ro') ;
plot(tms, +anglesBodyFrame(:, PHIL),'bo') ;
plot(tms, smooth_anglesMat_R(1,:),'r-','linewidth',1) ;
plot(tms, smooth_anglesMat_L(1,:),'b-','linewidth',1) ;
xlabel('Time [ms]') ;
ylabel('Stroke angle [deg]') ;
grid on ; box on ;
axis([tms(1) tms(end) 0 200]) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);

%----------------

ax2 = subplot(3,1,2) ;
hold on ;
plot(tms, anglesBodyFrame(:, THETAR),'ro') ;
plot(tms, anglesBodyFrame(:, THETAL),'bo') ;
plot(tms, smooth_anglesMat_R(2,:),'r-','linewidth',1) ;
plot(tms, smooth_anglesMat_L(2,:),'b-','linewidth',1) ;
xlabel('Time [ms]') ;
ylabel('Deviation angle [deg]') ;
grid on ; box on ;
axis([tms(1) tms(end) -10 60]) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);

ax3 = subplot(3,1,3) ;
hold on ;
plot(tms, anglesBodyFrame(:, ETAR),'ro') ;
plot(tms, anglesBodyFrame(:, ETAL),'bo') ;
plot(tms, smooth_anglesMat_R(3,:),'r-','linewidth',1) ;
plot(tms, smooth_anglesMat_L(3,:),'b-','linewidth',1) ;
xlabel('Time [ms]') ;
ylabel('Wing pitch angle [deg]') ;
grid on ; box on ;
axis([tms(1) tms(end) 0 200]) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);

linkaxes([ax1,ax2,ax3],'x')

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
set([hrl h_wa hpitch],'paperpositionmode','auto') ;
set(gca,'xlim',[tms(1) tms(end)]) ;
set(gca,'fontsize',14) ;
title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
print(gcf,[datapath 'phi_lr'],'-dpng','-r300') ;
set(gca,'xlim',[correctionTime(1) correctionTime(2)]) ;
print(gcf,[datapath 'phi_lr_zoom'],'-dpng','-r300') ;
set(gca,'xlim',[tms(1) tms(end)]) ;
saveas(gcf,[datapath 'phi_lr'],'fig') ;


figure(h_wa) ;sb1= subplot(3,1,1) ; sb2=subplot(3,1,2) ; sb3 = subplot(3,1,3) ;
sbvec = [sb1 sb2 sb3];
set(sbvec,'xlim',[tms(1) tms(end)]) ;

set(sbvec,'fontsize',14) ;
axes(sb1) ; title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;

print(gcf,[datapath 'wing_angles'],'-dpng','-r300') ;
set(sbvec,'xlim',[correctionTime(1) correctionTime(2)]) ;
print(gcf,[datapath 'wing_angles_zoom'],'-dpng','-r300') ;
set(sbvec,'xlim',[tms(1) tms(end)]) ;
saveas(gcf,[datapath 'wing_angles'],'fig') ;

