function [sp_y, smooth_y, rms_y, EstErr] = ...
    varySpline(x,y,EstErr_low,EstErr_high,N,maxRMS,plotFlag) 
%--------------------------------------------------------------------------
%This basically scans the space of splines, allows you to fix an upper
%bound on the rms error between fit and data, and then finds the second
%spline within this range that has minimal second derivative. This is sort
%of unnecessary right now (since you shouldn't scan regions above your max RMS
% err), but I would like to make it more sophisticated, so that it
% minimizes some weighted sum of 2nd derivative deviation and rms error. 
%
%Inputs:
%   x = input x data
%   y = input y data
%   EstErr_low = lower bound on spline error
%   EstErr_high = upper bound on spline error
%   N = number of points to sample in region [EstErr_low, EstErr_high]
%   maxRMS = upper bound on RMS error
%   plotFlag = should it show plots
%
%Outputs:
%   sp_y = spline 
%   smooth_y = smoothed data
%   rms_y = rms error
%   EstErr = input to mySplineSmooth that was used for the result
%
%----------------------------------------------------------------------------
EstErr_array = linspace(EstErr_low,EstErr_high,N) ;
plotColor = linspace(0,1,N) ;
rms_array = zeros(1,N) ;
meanSecondDeriv = zeros(1,N) ;
if plotFlag
    figure ;
end

for i = 1:N
    [sp_y_temp, smooth_y_temp, ~] = ...
        mySplineSmooth(x, y, EstErr_array(i)) ;
    
    rms_array(i) = sqrt(mean((smooth_y_temp - y).^2)) ;
    meanSecondDeriv(i) = mean(abs(fnval(fnder(sp_y_temp,2),x))) ;
    
    if plotFlag 
        s1 = subplot(1,4,1);
            hold on
            plot3(x, y, EstErr_array(i)*zeros(size(x)),'k.')
            plot3(x, smooth_y_temp, EstErr_array(i)*ones(size(x)),...
                'Color',[0 plotColor(i) 0], 'LineWidth',1)
            set(gca,'xlim',[min(x), max(x)])
            view([0 0 1])
            title('Spline and Data')
            xlabel('X')
            ylabel('Y')
        s2 = subplot(1,4,2) ;
            hold on
            plot3(x, fnval(fnder(sp_y_temp,1),x), EstErr_array(i)*ones(size(x)),...
                'Color',[0 plotColor(i) 0], 'LineWidth',2)
            set(gca,'xlim',[min(x), max(x)])
            view([0 0 1])
            title('First Derivative')
        s3 = subplot(1,4,3) ;
            hold on
            plot(EstErr_array(i), rms_array(i),'ko','MarkerSize',10,'MarkerFaceColor',...
                [0 plotColor(i) 0])
            title('RMS Error')
        s4 = subplot(1,4,4) ;
            hold on
            plot(EstErr_array(i), meanSecondDeriv(i),'ko','MarkerSize',10,'MarkerFaceColor',...
            [0 plotColor(i) 0])
            title('Mean Abs 2nd Derivative')
    end
end

temp_ind = find(rms_array <= maxRMS) ;
[minSecondDeriv, ~] = min(meanSecondDeriv(temp_ind)) ;
ind = find(meanSecondDeriv == minSecondDeriv) ;

EstErr = EstErr_array(ind) ;
rms_y = rms_array(ind) ;  
[sp_y, smooth_y, ~] = mySplineSmooth(x,y,EstErr) ; 

if plotFlag
    subplot(s1)
        plot3(x, smooth_y, EstErr_array(i)*ones(size(x)),...
                'Color',[.9 0 0], 'LineWidth',1.5)
    subplot(s2)
        plot3(x, fnval(fnder(sp_y,1),x), EstErr_array(i)*ones(size(x)),...
                'Color',[.9 0 0], 'LineWidth',2)
    subplot(s3)
        plot(EstErr, rms_y,'ko','MarkerSize',10,'MarkerFaceColor',[.9 0 0])
    subplot(s4)
        plot(EstErr, minSecondDeriv,'ko','MarkerSize',10,'MarkerFaceColor',[.9 0 0])
end

%disp('end'); 





