cd('D:\Janelia Flies\Analysis\Pitch Up\Expr_13_mov_002')
load('Expr13mov002_Data_manually_corrected.mat')

xlim = data.manualCorrRangeMS ;

pitchYlim = [10 70] ;
rollYlim = [-35 60] ;
yawYlim = [-15 45] ;

strokeYlim = [25 190] ;
rotYlim = [0 180] ;
devYlim = [-15 55] ;

numberSize = 12 ;
fontSize = 10 ;
plotWidth = 420 ;
plotHeight = 140 ;

%-----------------------------------------------------
hpitch =  bodyAnglePlot(data, true, false , false, false) ;
set(gca, 'xlim', xlim)
set(gca, 'ylim', pitchYlim)
set(gca, 'XTickLabel', []);
set(gca,'fontsize',10) ;
%set(gcf, 'Units', 'inches');
set(gcf, 'Position', [500 500 200 100]);
set(gcf, 'Color', 'none')
set(gcf, 'Renderer','opengl')

hyaw =  bodyAnglePlot(data, false, true , false, false) ;
set(gca, 'xlim', xlim)
set(gca, 'ylim', yawYlim)
set(gca, 'XTickLabel', []);
set(gca,'fontsize',10)
%set(gcf, 'Units', 'inches');
set(gcf, 'Position', [500 500 200 100]);
set(gcf, 'Color', 'none')
set(gcf, 'Renderer','opengl')

hroll =  bodyAnglePlot(data, false, false , true, false) ;
set(gca, 'xlim', xlim)
set(gca, 'ylim', rollYlim)
set(gca, 'XTickLabel', []);
set(gca,'fontsize',10)
%set(gcf, 'Units', 'inches');
set(gcf, 'Position', [500 500 200 100]);
set(gcf, 'Color', 'none')
set(gcf, 'Renderer','opengl')

hstroke = wingAnglePlot(data,0,0,0,strokeYlim) ; %stroke angle %[10 180]for 7.9
set(gca, 'XTickLabel', []);
%set(gca, 'YTick', [10:40:180]);
set(gca,'fontsize',numberSize)
%set(gcf, 'Units', 'inches');
set(gcf, 'Position', [500 500 plotWidth plotHeight]);
set(gcf, 'Color', 'w')
set(gcf, 'Renderer','opengl')
xlabel('')
title('')
ylabel('')
%ylabel({'Wing Stroke Angle'; '[deg]'},'fontsize',fontSize) ;
%ylabh = get(gca,'XLabel');
%set(ylabh,'Position',get(ylabh,'Position') + [0 1 0])
%xlabel({'Time [ms]'},'fontsize',.1) ;

hrot = wingAnglePlot(data,1,0,0,rotYlim) ;%wing pitch  %[15 180] for 7.9
set(gca,'fontsize',numberSize)
set(gca, 'XTickLabel', []);
%set(gca, 'YTick', [20:40:180]);
%set(gcf, 'Units', 'inches');
set(gcf, 'Position', [500 500 plotWidth plotHeight]);
set(gcf, 'Color', 'w')
set(gcf, 'Renderer','opengl')
set(gca,'ylim',[15 180])
xlabel('')
ylabel('')
%ylabel({'Wing Pitch Angle [deg]'},'fontsize',fontSize) ;
title('')

hdev = wingAnglePlot(data,2,0,0,devYlim) ; %wing deviation %[-5 40] for 7.9
set(gca, 'XTickLabel', []);
set(gca,'fontsize',numberSize)
%set(gcf, 'Units', 'inches');
set(gcf, 'Position', [500 500 plotWidth plotHeight]);
set(gcf, 'Color', 'w')
set(gcf, 'Renderer','opengl')
xlabel('')
title('')
ylabel('')
%ylabel({'Wing Deviation Angle [deg]'},'fontsize',fontSize) ;
xlabel('') ;

%{
plot2svg('bodypitchangle.svg',hpitch,'png')
plot2svg('bodyyawangle.svg',hyaw,'png')
plot2svg('bodyrollangle.svg',hroll,'png')
plot2svg('strokeangle.svg',hstroke,'png')
plot2svg('wingpitchangle.svg',hrot,'png')
plot2svg('deviationangle.svg',hdev,'png')
%}