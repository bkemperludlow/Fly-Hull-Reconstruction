function [i1, i2] = chooseCutWindow(angleR,angleL,t,window,backFlipTimes,fwdFlipTimes)

patchColor = [1 1 1]*.8 ;

figure ;
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [2 4 12 6]);
hold on
plot(t*1000, angleR, 'ro-')
plot(t*1000, angleL, 'bo-')
plot(t(window(1))*1000, angleR(window(1)), 'ksq','MarkerSize',10,'LineWidth',2)
plot(t(window(2))*1000, angleR(window(2)), 'ksq','MarkerSize',10,'LineWidth',2)
plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
set(gca,'xlim',[t(window(1)-20)*1000, t(window(2)+20)*1000]) ;

try
    [x1, y1] = ginput(1) ;
    [x2, y2] = ginput(1) ;
catch 
    disp('Need to click valid graph points')
    return ;
end

if x1 > x2
    x_temp = x2 ;
    y_temp = y2 ;
    x2 = x1 ;
    y2 = y1 ;
    x1 = x_temp ;
    y1 = y_temp ;
end

diff1 = abs(t - x1/1000) ;
[~, i1] = min(diff1) ;
diff2 = abs(t - x2/1000) ;
[~, i2] = min(diff2) ;



