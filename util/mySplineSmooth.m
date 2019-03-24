function  [sp, smoothed_y, err]  = mySplineSmooth(x, y, desiredErr, weights)
%
% [sp, smoothed_y, err]  = mySplineSmooth(x, y, desiredErr )
%
% spline smooth y(x), such the the mean-square error of the smoothed data
% is desiredErr. 
%
% input parameters:
%    x and y     : two vectors that define the function y(x). 
%                  x does not have to be evenly spaced.
%    desiredErr  : the desired mean-square error of the fit (scalar).
%    weights     : gives weighting for different points (size of x)
%                  see spaps documentation for details                   
%
% output parameters:
%    sp         : the matlab spline object given by spaps
%    smoothed_y : the value of the spline calculated for x.
%    err        : the resulting mean-square error, decined as:
%                 sqrt(mean((y-smoothed_y).^2))
%
% Example:
% 
% N=200 ; 
% noiseAmp = 0.1 ;
% x = linspace(0,2*pi,N) ; 
% y = sin(x) + (rand(1,N)-0.5)*noiseAmp ;
% [sp, smoothed_y, err]  = mySplineSmooth(x, y, noiseAmp/3 ) ;
% figure, plot(x,y,'ko','markerfacecolor','c') ; hold on ;
% plot(x,smoothed_y,'r-','linewidth',2) ;
% set(gca,'xlim',[x(1) x(end)]); box on ; grid on ; xlabel('x') ; ylabel('y') ;

N = length(x) ;
if (length(y)~=N)
    error('x and y shoud have the same length') ;
end

% reshape y to have the same dimension as x
if (size(x)~=size(y))
    y=y';
end

dx  = mean(diff(x)) ;
tol = N * desiredErr^2 * dx ;

if exist('weights','var')
    [sp, smoothed_y] =  spaps(x, y, tol,weights) ;
else    
    [sp, smoothed_y] =  spaps(x, y, tol) ;
end

% verify fit mean-square error
%err = ( mean( (y-smoothed_y).^2))^0.5 ;
err = [];

end

