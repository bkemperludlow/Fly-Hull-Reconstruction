%--------------------------------------------------------------------------
% main function to calculate body and wing Euler angles for a single fly
% video. calculates angles in both lab and body frame.
%--------------------------------------------------------------------------
function data = calcAnglesMain(dataPath, largePertFlag, saveFlag, plotFlag,...
    suffixStr)
%--------------------------------------------------------------------------
%% params and inputs
if ~exist('dataPath','var') || isempty(dataPath)
    %mFilePath = mfilename('fullpath') ;
    %[filterSpec,~,~] = fileparts(mFilePath) ;
    %filterSpec =  'D:\Fly Data\VNC MN Chrimson\03_06112017\Analysis\Unsorted\Expr_3_mov_000\' ; 
    filterSpec = 'D:\Fly Data\VNC MN Chrimson\21_03102018\Analysis\Unsorted\Expr_21_mov_037\' ;
    %filterSpec = 'D:\Fly Data\VNC MN Chrimson\52_09112019\Analysis\Unsorted\Expr_52_mov_008\' ; 
    filepath_cell = uipickfiles('FilterSpec',filterSpec,...
        'Type',{'*.mat', 'MAT-files' },'NumFiles',1) ;
    try
        dataPath = filepath_cell{1} ;
    catch
        disp('No data file selected')
        data = [] ; 
        return
    end
    
end
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false ;
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = true ; % true
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true ;
end
if ~exist('suffixStr','var') || isempty(suffixStr)
   suffixStr = '' ;  
end

defineConstantsScript
%smoothingParams = setSmoothingParams ;
pulseStartMS = 0 ;
%pulseLengthMS = 7 ;
%dt = 1/data.params.fps ; 
patchColor = [1 1 1 ] * 0.8;

%--------------------------------------------------------------------------
%% load in data file
[datapath, datafilename, ~] = fileparts(dataPath) ; 

if contains(datafilename,'Expr_')
    % this is the case for non-manually corrected data
    datafilename_split = strsplit(datafilename,'_') ;
    exprNum = str2double(datafilename_split{2});
    movNum = str2double(datafilename_split{4}) ;
    sufStr = datafilename_split{end} ; 
    if strcmp(suffixStr,'') && ~(strcmp(sufStr,'test') || strcmp(sufStr,'results'))
        suffixStr = ['_' sufStr] ; 
    end
else
    try
        % this is the case for manually corrected data
        expr_idx = strfind(datafilename,'Expr') ;
        mov_idx = strfind(datafilename,'mov') ;
        exprNum = str2double(datafilename((expr_idx+4):(mov_idx-1))) ;
        movNum = str2double(datafilename((mov_idx+3):(mov_idx+5))) ;
    catch
        % if filename is nonstandard (doesn't really matter)
        exprNum = nan ; 
        movNum = nan ; 
    end
end

% load data
data = importdata(dataPath) ;

%--------------------------------------------------------------------------
%% read in data fields regarding manual correction (if they exist)
if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end
% * new as of july 2019: we try to identify frames in which wings are
% severly clipped, in which case they're probably best to ignore
% if (isfield(data,'clippedWingFlag'))
%    ignoreFrames = sort([ignoreFrames, find(data.clippedWingFlag)']) ; 
% end
% manual correction time, if mc has occured. gives plot bounds
if (isfield(data,'correctionTime'))
    correctionTime=data.correctionTime;
elseif(isfield(data,'manualCorrRangeMS'))
    correctionTime=data.manualCorrRangeMS;
else
    correctionTime = [-10, 50] ;
end

%--------------------------------------------------------------------------
%% calculate raw angles
[anglesLabFrame, anglesBodyFrame, t, ~, ~, ~, smoothed_rho, rho_t, ...
    rho_samp, rotM_YP, rotM_roll, largePertFlag ] = ...
        calcAnglesRaw_Sam(data, plotFlag,largePertFlag) ;

%--------------------------------------------------------------------------
%% unwrap and spline smooth wing angles
%------------------------------------------------
% right stroke angle
phiR = -anglesBodyFrame(:, PHIR) ;
ignoreIndR = unique([find(isnan(phiR))' ignoreFrames]) ;
%phiR = phiR + 360 ;

for i = 1:length(phiR)
    while phiR(i) < -90
        phiR(i) = phiR(i) + 360 ;
    end
    while phiR(i) > 270
        phiR(i) = phiR(i) - 360 ;
    end
end

if (~isempty(ignoreIndR))
    phiR(ignoreIndR) = NaN ;
end

%-----------------------------------------------
% left stroke angle
phiL = +anglesBodyFrame(:, PHIL) ;
ignoreIndL = unique([find(isnan(phiL))'  ignoreFrames])  ;

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

% hampel filter to remove outliers
[~, hampelR] = hampel(phiR, 7,2) ;
[~, hampelL] = hampel(phiL, 7,2) ;

phiR(hampelR) = NaN ;
phiL(hampelL) = NaN ;

%--------------------------------------------------------------------------
%% get wing flip times
[fwdFlipTimesR, backFlipTimesR, fwdFlipIndR, backFlipIndR, fwdFlipPhiR,...
    backFlipPhiR, badIndicesR] = findWingFlipTimes_mk3 (t, phiR, plotFlag);
[fwdFlipTimesL, backFlipTimesL, fwdFlipIndL, backFlipIndL, fwdFlipPhiL,...
    backFlipPhiL, badIndicesL] = findWingFlipTimes_mk3 (t, phiL, plotFlag);

phiR(badIndicesR) = NaN ;
phiL(badIndicesL) = NaN ;

%--------------------------------------------------------------------------
%% store data in structure
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

% 4/23/20: want to include additional info: rotation matrices, 
% largePertFlag, ExprNum, MovNum, smoothed angles
data.rotM_YP = rotM_YP ; 
data.rotM_roll = rotM_roll ; 
data.largePertFlag = largePertFlag ; 

data.ExprNum = exprNum ; 
data.MovNum = movNum ; 

%--------------------------------------------------------------------------
%% find stroke amplitude
phiRfwdFlip = fwdFlipPhiR' ;  % fnval(sp_phiR_low, fwdFlipTimesR) ;
phiRbckFlip = backFlipPhiR' ; % fnval(sp_phiR_low, backFlipTimesR) ;
phiLfwdFlip = fwdFlipPhiL' ;  % fnval(sp_phiL_low, fwdFlipTimesL) ;
phiLbckFlip = backFlipPhiL' ; % fnval(sp_phiL_low, backFlipTimesL) ;

% right wing
phir_comb = [phiRfwdFlip, phiRbckFlip] ;
phir_t    = [fwdFlipTimesR, backFlipTimesR] ;
rmat = [phir_t', phir_comb'] ;
rmat = sortrows(rmat,1) ;

phir_amp = abs(diff(rmat(:,2))) ;
phir_amp_t = rmat(1:end-1,1) + diff(rmat(:,1))/2 ;
mid_stroke_r = rmat(1:end-1,2) + diff(rmat(:,2))/2 ;

% left wing
phil_comb = [phiLfwdFlip, phiLbckFlip] ;
phil_t    = [fwdFlipTimesL, backFlipTimesL] ;
lmat = [phil_t', phil_comb'] ;
lmat = sortrows(lmat,1) ;

phil_amp = abs(diff(lmat(:,2))) ;
phil_amp_t = lmat(1:end-1,1) + diff(lmat(:,1))/2 ;
mid_stroke_l = lmat(1:end-1,2) + diff(lmat(:,2))/2 ;

% add stroke amplitude and associated times to struct
data.phil_amp_t = phil_amp_t ;
data.phir_amp_t = phir_amp_t ;
data.phiL_amp = phil_amp ;
data.phiR_amp = phir_amp ;

% -------------------------------------------------------------------------
%% smooth angles (body and wing)
% --------------------------------------------
% smooth wing angles (lab and body frames)
[~, smoothAnglesMatR_Lab, ~, ~, ~ ] = smoothWingAngles(data, 'R','Lab') ;
[~, smoothAnglesMatL_Lab, ~, ~, ~ ] = smoothWingAngles(data, 'L','Lab') ;
[~, smoothAnglesMatR_Body, ~, ~, ~ ] = smoothWingAngles(data, 'R','Body') ;
[~, smoothAnglesMatL_Body, ~, ~, ~ ] = smoothWingAngles(data, 'L','Body') ;
% make sure phiR is negative in body frame
if (mode(sign(smoothAnglesMatR_Body(1,:))) > 0)
    smoothAnglesMatR_Body(1,:) = -1.*smoothAnglesMatR_Body(1,:) ; 
end

% --------------------------------------------------------------------
% smooth body angles (just in lab frame -- not defined in body frame)
[pitchSmooth, yawSmooth, rollSmooth] = smoothBodyAngles(data,largePertFlag) ;

% ------------------------------------
% create arrays for smoothed angles
% NB: need to take transpose for wing angle mats
anglesLabFrameSmooth = [yawSmooth, pitchSmooth, smoothAnglesMatR_Lab', ...
    smoothAnglesMatL_Lab', rollSmooth] ; 
anglesBodyFrameSmooth = zeros(data.Nimages, 8);
anglesBodyFrameSmooth(:,[PHIR, THETAR, ETAR, PHIL, THETAL, ETAL]) = ...
    [smoothAnglesMatR_Body', smoothAnglesMatL_Body'] ; 
    
% --------------------------------------
% add to data struct
data.anglesLabFrameSmooth = anglesLabFrameSmooth ; 
data.anglesBodyFrameSmooth = anglesBodyFrameSmooth ; 

%--------------------------------------------------------------------------
%% save results?
if saveFlag
    save(dataPath,'data')
end

%--------------------------------------------------------------------------
%% plot results?
if plotFlag
    % smooth wing angles
    tms = t * 1000 ;
    
    %----------------------------------------------------------------------
    % R/L stroke angle with body pitch
    hrl = figure('position',[140 550 1100 420])  ; hold on ;
    
    plot(tms, -anglesBodyFrame(:, PHIR),'ro') ;
    plot(tms, +anglesBodyFrame(:, PHIL),'bo') ;
    plot(phir_amp_t*1000, mid_stroke_r,'rs-','linewidth',2) ;
    plot(phil_amp_t*1000, mid_stroke_l,'bs-','linewidth',2) ;
    plot(tms, anglesLabFrame(:,BETA),'.','color',[0.2 0.8 0.2]*0) ;
    plot(tms, anglesLabFrameSmooth(:,BETA),'-','color',[0.2 0.8 0.2]*0) ;
    plot(tms, -anglesBodyFrameSmooth(:, PHIR),'r-') ;
    plot(tms, +anglesBodyFrameSmooth(:, PHIL),'b-') ;
    xlabel('Time [ms]') ;
    ylabel('Stroke angle [deg]') ;
    grid on ; box on ;
    axis([tms(1) tms(end) 0 200]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    
    
    plot(fwdFlipTimesR*1000, phiRfwdFlip,'k^','markerfacecolor','r') ;
    plot(backFlipTimesR*1000, phiRbckFlip,'kv','markerfacecolor','r') ;
    
    plot(fwdFlipTimesL*1000, phiLfwdFlip,'k^','markerfacecolor','b') ;
    plot(backFlipTimesL*1000, phiLbckFlip,'kv','markerfacecolor','b') ;
    
    %----------------------------------------------------------------------
    % wing angles
    h_wa = figure('position',[140 100 1100 800])  ;
    title('Wing Angles')
    
    %----------------
    % stroke
    ax1 = subplot(3,1,1) ;
    hold on ;
    plot(tms, -anglesBodyFrame(:, PHIR),'ro') ;
    plot(tms, +anglesBodyFrame(:, PHIL),'bo') ;
    plot(tms, -anglesBodyFrameSmooth(:, PHIR),'r-','linewidth',1) ;
    plot(tms, +anglesBodyFrameSmooth(:, PHIL),'b-','linewidth',1) ;
    xlabel('Time [ms]') ;
    ylabel('Stroke angle [deg]') ;
    grid on ; box on ;
    axis([tms(1) tms(end) 0 200]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    
    %----------------
    % deviation
    ax2 = subplot(3,1,2) ;
    hold on ;
    plot(tms, anglesBodyFrame(:, THETAR),'ro') ;
    plot(tms, anglesBodyFrame(:, THETAL),'bo') ;
    plot(tms, anglesBodyFrameSmooth(:, THETAR),'r-','linewidth',1) ;
    plot(tms, anglesBodyFrameSmooth(:, THETAL),'b-','linewidth',1) ;
    xlabel('Time [ms]') ;
    ylabel('Deviation angle [deg]') ;
    grid on ; box on ;
    axis([tms(1) tms(end) -10 60]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    
    %----------------
    % pitch
    ax3 = subplot(3,1,3) ;
    hold on ;
    plot(tms, anglesBodyFrame(:, ETAR),'ro') ;
    plot(tms, anglesBodyFrame(:, ETAL),'bo') ;
    plot(tms, anglesBodyFrameSmooth(:, ETAR),'r-','linewidth',1) ;
    plot(tms, anglesBodyFrameSmooth(:, ETAL),'b-','linewidth',1) ;
    xlabel('Time [ms]') ;
    ylabel('Wing pitch angle [deg]') ;
    grid on ; box on ;
    axis([tms(1) tms(end) 0 200]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS,...
        fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    
    linkaxes([ax1,ax2,ax3],'x')
    
    %----------------------------------------------------------------------
    % raw pitch angle
    hpitch = figure ;
    plot(tms, anglesLabFrameSmooth(:,BETA),'-','color',0*[0.2 0.8 0.2]) ;
    plot(tms, anglesLabFrame(:,BETA),'o','color',[0.2 0.8 0.2]) ;
    xlabel('Time [ms]') ;
    ylabel('Body pitch angle [deg]') ;
    grid on ; box on ;
    set(gca,'xlim',[tms(1) tms(end)]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, ...
        fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    
    %----------------------------------------------------------------------
    % roll angle
    hroll = figure ;
    hold on ;
    plot(rho_t*1000, rho_samp,'ko','markerfacecolor','g') ;
    plot(tms, smoothed_rho,'k-') ;
    plot(tms, anglesLabFrameSmooth(:,RHO),'k--')
    set(gca,'fontsize',14) ;
    title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum)]) ;
    xlabel('Time [ms]') ;
    ylabel('Body roll angle [deg]') ;
    grid on ; box on ;
    set(gca,'xlim',[tms(1) tms(end)]) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000-pulseStartMS, fwdFlipTimesR*1000-pulseStartMS, patchColor, true);
    if saveFlag
        savefig(gcf,fullfile(datapath, ['roll_angle' suffixStr '.fig'])) ;
        print(gcf,fullfile(datapath, ['roll_angle' suffixStr '.png']),...
            '-dpng','-r300') ;
    end
    
    
    %axis([tms(1) tms(end) 30 70]) ;
    
    %% save plots?
    if saveFlag
        figure(hrl) ;
        set([hrl h_wa hpitch],'paperpositionmode','auto') ;
        set(gca,'xlim',[tms(1) tms(end)]) ;
        set(gca,'fontsize',14) ;
        title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum,'%03d')]) ;
        print(gcf,fullfile(datapath, ['phi_lr' suffixStr '.png']),...
            '-dpng','-r300') ;
        set(gca,'xlim',[correctionTime(1) correctionTime(2)]) ;
        print(gcf,fullfile(datapath, ['phi_lr_zoom' suffixStr '.png']),...
            '-dpng','-r300') ;
        set(gca,'xlim',[tms(1) tms(end)]) ;
        savefig(gcf,fullfile(datapath, ['phi_lr' suffixStr '.fig'])) ;
        
        
        figure(h_wa) ;
        sb1= subplot(3,1,1) ; 
        sb2=subplot(3,1,2) ; 
        sb3 = subplot(3,1,3) ;
        sbvec = [sb1 sb2 sb3];
        set(sbvec,'xlim',[tms(1) tms(end)]) ;
        
        set(sbvec,'fontsize',14) ;
        axes(sb1) ;
        title(['Expr ' num2str(exprNum) ' Movie ' num2str(movNum,'%03d')]) ;

        print(gcf,fullfile(datapath, ['wing_angles' suffixStr '.png']),...
            '-dpng','-r300') ;
        set(sbvec,'xlim',[correctionTime(1) correctionTime(2)]) ;
        print(gcf,fullfile(datapath, ['wing_angles_zoom' suffixStr '.png']),...
            '-dpng','-r300') ;
        set(sbvec,'xlim',[tms(1) tms(end)]) ;
        savefig(gcf,fullfile(datapath, ['wing_angles' suffixStr '.fig'])) ;
    end
end
end