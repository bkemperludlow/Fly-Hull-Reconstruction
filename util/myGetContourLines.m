% -------------------------------------------------------------------------
% function to get x, y coordinates of contour lines from contour matrix
% (output of either contourc or contour plotting functions).
%
% Heavily borrowed from:
% https://www.mathworks.com/matlabcentral/fileexchange/74010-getcontourlinecoordinates
% -------------------------------------------------------------------------
function [vertex_list, level_list, num_list] = myGetContourLines(cm)
% Determine if input is handle or matrix; get matrix.
if ishandle(cm)
    cm = cm.ContourMatrix;
end

% Set up while loop
cmSize = size(cm,2);   	% number of columns in ContourMatrix
cmWindow = [0,0];      	% [start,end] index of moving window
contourArray = {};     	% Store the (x,y) coordinates of each contour line
% Extract coordinates of each contour line
while cmWindow(2) < cmSize
    cmWindow(1) = cmWindow(2) + 1;
    cmWindow(2) = cmWindow(2) + cm(2,cmWindow(1)) + 1;
    contourArray(end+1) = {cm(:,cmWindow(1):cmWindow(2)).'};  %#ok<AGROW>
end

% Separate the level, count, and coordinates.
level_list = cellfun(@(c)c(1,1),contourArray).';
%numCoord = cellfun(@(c)c(1,2),contourArray).';
vertex_list = cellfun(@(c)c(2:end,:),contourArray,'UniformOutput',false);

% Sort by level (just in case Matlab doesn't)
[level_list, sortIdx] = sort(level_list);

vertex_list = vertex_list(sortIdx) ; 
num_list = 1:length(vertex_list) ; 


end