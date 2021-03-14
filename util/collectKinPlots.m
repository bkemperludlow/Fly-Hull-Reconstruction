% -------------------------------------------------------------------------
% function to find plots of specific kinematic variables from movies in a 
% single experiment and copy them all to a single folder for easy viewing
%
% INPUTS:
% --------
%   rootPath: string path to directory (e.g. '...\75_14062019\']
%   plotType: string that determines which types of plots to copy over. 
%       can be any of 'wing_angles', 'roll_angle', or 'phi_lr'
% -------------------------------------------------------------------------
function [] = collectKinPlots(rootPath, plotType, fileExt, excludeFlag)
% ----------------
% inputs/params:
if ~exist('fileExt','var') || isempty(fileExt)
    [~,~,fileExt] = fileparts(plotType) ; 
    if isempty(fileExt)
        fileExt = '.png' ; 
    end
end
if ~exist('excludeFlag','var') || isempty(excludeFlag)
    excludeFlag = true ; 
end
pathStruct = generatePathStruct(rootPath) ; 
% -------------------------------------------------------------
% get directory for all instances of plot type within folder
[~, plotTypeName, ~] = fileparts(plotType) ; 
dirPath = fullfile(rootPath,'**','*',[plotTypeName fileExt]) ;
plotDir = dir(dirPath) ;

% exclude plots that are separated from main data set
if excludeFlag
   exclude_idx = arrayfun(@(x) contains(x.folder, 'Exclude',...
       'IgnoreCase',0), plotDir) ;
   plotDir = plotDir(~exclude_idx) ; 
end
N_plots = length(plotDir) ; 

if (N_plots < 1)
   fprintf('Could not find %s plots in folder %s of format %s \n',...
       plotTypeName, rootPath, fileExt)
   return
end

% -----------------------------------------------------------
% if there are plots, create folder to move them to
newPlotFolder = fullfile(pathStruct.root, [plotTypeName '_plots']) ;
if exist(newPlotFolder,'dir') < 1
    mkdir(newPlotFolder)
end

% -------------------------------------------------------------------------
% get experiment and movie numbers so that when we copy plots we can rename
% them (by default, the kinematics plots don't include experiment and movie
% identifiers)
ExprNum = pathStruct.ExprNum ; 
MovNumStrList = arrayfun(@(x) x.folder(end-2:end), plotDir,...
    'UniformOutput', 0) ;

% check that we're actually getting a movie number
good_idx = cellfun(@(y) ~isnan(str2double(y)), MovNumStrList) ; 
if any(~good_idx)
    plotDir = plotDir(good_idx) ; 
    MovNumStrList = arrayfun(@(x) x.folder(end-2:end), plotDir,...
        'UniformOutput', 0) ;
    N_plots = length(plotDir) ; 
end

% ------------------------------
% copy over plots
for i = 1:N_plots
   prefixStr = ['Expr_' num2str(ExprNum) '_mov_' MovNumStrList{i} '_'] ;
   fname_old = fullfile(plotDir(i).folder, plotDir(i).name) ; 
   fname_new = fullfile(newPlotFolder, strcat(prefixStr, plotTypeName, fileExt)) ; 
   
   [status, msg] = copyfile(fname_old, fname_new) ; 
   if status
       fprintf('Successfully copied %s \n to %s \n', fname_old, fname_new)
   else
       fprintf('Error: %s \n', str(msg)) 
   end
end

end