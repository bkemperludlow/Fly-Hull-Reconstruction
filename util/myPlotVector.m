function h = myPlotVector(v1,v2,color, width, linestyle)
% plots the vector v2 such that its base point is v1
% using the given color 
if (~exist('width','var'))
    width = 2 ;
end
if (~exist('linestyle','var'))
    linestyle = '-';
end

v3 = v1+v2 ;
h = plot3( [v1(1) v3(1)], [v1(2) v3(2)], [v1(3) v3(3)],'Color',color,...
    'linewidth',width, 'LineStyle', linestyle) ;
end