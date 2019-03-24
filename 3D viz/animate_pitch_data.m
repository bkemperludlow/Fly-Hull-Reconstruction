dataPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis\' ;
savePath = 'F:\Dropbox\Paper Manuscripts\b1 paper\figure_3\' ;
% dataPath = 'H:\Fly Data\Opto Silencing\01_26022018\Analysis\Unsorted\Expr_1_mov_012\' ;
% savePath = 'F:\Dropbox\Paper Manuscripts\b1 paper\figure_5\' ;

b1_talk_plot_pref

ExprNum = 5 ;
MovNum = 2 ;

camView = 'XZ' ;
flipImFlag = true ;
saveFinalFlag = true ;
FPS = 30 ;
% only necessary if you want images:

analysisPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Up\' ;
cinRoot = 'G:\Janelia Flies\kir2.1 flies round 2' ;

figPosition = [7.0729    5.8021    7.1771    2.0313] ;
%imFigPosition =  [7.0729    5.8021    3.25    2.5] ;

lineWidth_thick = 2.5 ;
lineWidth_thin = 1.5 ;

color_cell = {ctrl_kir_pitch_color, mb204b_ctrl_pitch_color,...
    mb204b_kir_pitch_color} ;


if ismember(ExprNum,[3,5,6,14,15])
    plotColor = color_cell{1} ;
elseif ismember(ExprNum,[2,7,10])
    plotColor = color_cell{2} ;
elseif ismember(ExprNum,[1,8,11])
    plotColor = color_cell{3} ;
else
    plotColor = [0 0 0] ;
end

%plotColor = [117,107,177]/255 ;
%--------------------------------------------------------------------------

controller_fit_struct_all = importdata(fullfile(dataPath,'pitch_controller_fits_all.mat')) ;
%controller_fit_struct_all = importdata(fullfile(dataPath,'controller_fit_struct.mat')) ;

data_ind = ([controller_fit_struct_all.ExprNum] == ExprNum) & ...
    ([controller_fit_struct_all.MovNum] == MovNum) ;

if sum(data_ind) ~= 1
    keyboard
end

controller_curr = controller_fit_struct_all(data_ind) ;

fwdFlipTimes = controller_curr.fwdFlipTimes ;
c_pitch = controller_curr.c_pitch ;
deltaPhiFront = controller_curr.deltaPhiFront ;

theta_0 = c_pitch(0) ;
t = fwdFlipTimes(1): 0.000125 : fwdFlipTimes(end) ;

deltaBodyPitch_cont = c_pitch(t) - theta_0 ; %continuous

im_struct = get_fly_snapshots(cinRoot, ExprNum, MovNum, camView,...
    analysisPath, t,savePath, false, false)  ;

%--------------------------------------------------------------------------
%====================
% phi front figure
%====================
h_phiFront = figure('PaperPositionMode','auto','units',figUnits,...
    'Position',figPosition,'Color','w') ;
hold on

ax_phiFront = gca ;
ylim_phiFront = [min(deltaPhiFront)-5 , max(deltaPhiFront)+5] ;
set(ax_phiFront, 'LineWidth', axisLineWidth)
set(ax_phiFront, 'xlim', 1000*[t(1), t(end)])
set(ax_phiFront,'ylim',ylim_phiFront)
set(ax_phiFront, 'fontname',fontName,'fontsize',axisFontSize)
xlabel('Time [ms]', 'fontname',fontName,'fontsize',labelFontSize)
ylabel('\Delta \Phi_{front} [deg]', 'fontname',fontName,'fontsize',labelFontSize)

saveName_phiFront = [savePath 'Expr_' num2str(ExprNum) '_mov_' ...
    num2str(MovNum,'%03d') '_phiFront_animated.mp4'] ;
vidWriter_phiFront = VideoWriter(saveName_phiFront,'MPEG-4') ;
vidWriter_phiFront.FrameRate = FPS ;

%====================
% body pitch figure
%====================
h_pitch = figure('PaperPositionMode','auto','units',figUnits,...
    'Position',figPosition,'Color','w') ;
hold on

ax_pitch = gca ;
ylim_pitch = [min(deltaBodyPitch_cont)-5 , max(deltaBodyPitch_cont)+5] ;
set(ax_pitch, 'LineWidth', axisLineWidth)
set(ax_pitch, 'xlim', 1000*[t(1), t(end)])
set(ax_pitch,'ylim',ylim_pitch)
set(ax_pitch, 'fontname',fontName,'fontsize',axisFontSize)
xlabel('Time [ms]', 'fontname',fontName,'fontsize',labelFontSize)
ylabel('\Delta Body Pitch [deg]', 'fontname',fontName,'fontsize',labelFontSize)

saveName_pitch = [savePath 'Expr_' num2str(ExprNum) '_mov_' ...
    num2str(MovNum,'%03d') '_pitch_animated.mp4'] ;
vidWriter_pitch = VideoWriter(saveName_pitch,'MPEG-4') ;
vidWriter_pitch.FrameRate = FPS ;
%====================
% image figure
%====================
h_images = figure('PaperPositionMode','auto','units',figUnits,'Color','w') ;
ax_im = gca ;
saveName_im = [savePath 'Expr_' num2str(ExprNum) '_mov_' ...
    num2str(MovNum,'%03d') '_images.mp4'] ;
vidWriter_im = VideoWriter(saveName_im,'MPEG-4') ;
vidWriter_im.FrameRate = FPS ;

%--------------------------------------------------------------------------
open(vidWriter_phiFront) ;
open(vidWriter_pitch) ;
open(vidWriter_im) ;

for t_ind = 1:length(t)
    
    curr_t = t(t_ind) ;
    if curr_t >= 0
        
        pulse_DT = min([1000*curr_t, pulseDurMagnetic]) ;
        tsfvec = [0 pulse_DT pulse_DT 0 0] ;
        %patchColor = [1 1 1 ] * 0.8;
        
        avec_phiFront = [ylim_phiFront(1) ylim_phiFront(1) ylim_phiFront(2) ...
            ylim_phiFront(2) ylim_phiFront(1)] ;
        hf_phiFront = fill(tsfvec , avec_phiFront,'y','Parent',ax_phiFront) ;
        set(hf_phiFront,'facecolor',magnetic_pulse_color,'facealpha',pulseAlpha,'edgecolor','none') ;
        set(hf_phiFront,'HandleVisibility','off')
        
        avec_pitch = [ylim_pitch(1) ylim_pitch(1) ylim_pitch(2) ...
            ylim_pitch(2) ylim_pitch(1)] ;
        hf_pitch = fill(tsfvec , avec_pitch,'y','Parent',ax_pitch) ;
        set(hf_pitch,'facecolor',magnetic_pulse_color,'facealpha',pulseAlpha,'edgecolor','none') ;
        set(hf_pitch,'HandleVisibility','off')
        
    end
    
    hh_pitch = plot(1000*t(1:t_ind),deltaBodyPitch_cont(1:t_ind),'Color',plotColor,...
        'LineWidth',lineWidth_thick,'Parent',ax_pitch) ;
    
    hh_data = errorbar(1000*fwdFlipTimes(fwdFlipTimes <= curr_t),...
        deltaPhiFront(fwdFlipTimes <= curr_t),...
        2*ones(size(deltaPhiFront(fwdFlipTimes <= curr_t))),...
        'ko','markerfacecolor','k','Parent',ax_phiFront) ;
    
    hh_im = imshow(im_struct(t_ind).image,'Parent',ax_im) ;
    
    % get frames
    frame_pitch = getframe(h_pitch) ;
    writeVideo(vidWriter_pitch,frame_pitch) ;
    
    frame_phiFront = getframe(h_phiFront) ;
    writeVideo(vidWriter_phiFront,frame_phiFront) ;
    
    frame_im = getframe(ax_im) ;
    writeVideo(vidWriter_im,frame_im) ;
    
end

close(vidWriter_phiFront) ;
close(vidWriter_pitch) ;
close(vidWriter_im) ;
%--------------------------------------------------------------------------
if saveFinalFlag
    K_i = controller_curr.K_i ;
    K_p = controller_curr.K_p ;
    deltaT = controller_curr.deltaT ;
    K = controller_curr.K ;
    
    
    theta_0 = c_pitch(0) ;
    t = linspace(fwdFlipTimes(1), fwdFlipTimes(end), 100) ;
    
    deltaBodyPitch_forController = c_pitch(t - deltaT) - theta_0 ; %continuous
    pitchVelocity_forController = differentiate(c_pitch, t - deltaT) ;
    controlPred = K_i * deltaBodyPitch_forController + K_p * pitchVelocity_forController + K ;
    
    [upper, lower] = get_controllerFit_CI_pitch( controller_curr) ;
    
    h_fill = fill([1000*t,fliplr(1000*t)],[lower',fliplr(upper')],...
        plotColor,'linestyle','none','Parent',ax_phiFront);
    h_fill.FaceAlpha = 0.2 ;
    
    h_pi = plot(1000*t,controlPred,'Color',plotColor,...
        'LineWidth',lineWidth_thick,'Parent',ax_phiFront) ;
    h_p = plot(1000*t,K_p * pitchVelocity_forController,'Color',p_term_color,...
        'LineWidth',lineWidth_thin,'LineStyle',p_term_lineStyle,...
        'Parent',ax_phiFront) ;
    h_i = plot(1000*t,K_i * deltaBodyPitch_forController,'Color',i_term_color,...
        'LineWidth',lineWidth_thin,'LineStyle',i_term_lineStyle,...
        'Parent',ax_phiFront) ;
    
    saveNameFinal = [savePath 'Expr_' num2str(ExprNum) '_mov_' ...
        num2str(MovNum,'%03d') '_phiFront_animated_finalFrame'] ;
    
    savefig(h_phiFront, [saveNameFinal '.fig'])
    print(h_phiFront, [saveNameFinal '.svg'],'-dsvg')
    print(h_phiFront, [saveNameFinal '.png'],'-dpng','-r300')
    
end



