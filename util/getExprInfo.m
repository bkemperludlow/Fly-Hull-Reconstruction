% -------------------------------------------------------------------------
% function to grab experiment information (e.g. pulse timings, genotype,
% cine info) for a given experiment. assumes this info has been stored in
% the "*catalog.csv" for the experiment parent folder.
%
% if not, see the the "tempGetCineFrameRange.m" function in the
% "under_construction" folder to help with getting cine info -- the rest
% could likely be obtained from the README, but entry is non-uniform and if
% we keep updating manually as we go it should be okay
% -------------------------------------------------------------------------
function exprInfoStruct = getExprInfo(pathStruct, saveFlag, loadFlag)
% ----------------------------------------
%% inputs
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = true ; % save results?
end
if ~exist('loadFlag','var') || isempty(loadFlag)
    loadFlag = true ; % load results if they already exist?
end
% --------------------------------------
%% read out basic info from path struct
exprPath = pathStruct.root ;
ExprNum = pathStruct.ExprNum ;

% also initialize output struct
exprInfoStruct = struct() ;

% ------------------------------------------------------------------
%% see if exprInfoStruct already exists, in which case just load it
savePathFull = fullfile(exprPath, 'exprInfoStruct.mat') ; 
if exist(savePathFull,'file') && loadFlag
    % if we can find file, load and return
    exprInfoStruct = importdata(savePathFull) ; 
    return
end
% --------------------------------------------------------
%% try to locate catalog file for experiment parent folder
parentFolder = fullfile(exprPath, '..') ;
catalogDir = dir(fullfile(parentFolder,'*_catalog.csv')) ;

% check that we only get one
if length(catalogDir) ~= 1
    fprintf('Error: could not locate catalog file for: \n %s \n', exprPath)
    return
end

% if we have found catalog file, load it
catalog_fn = fullfile(catalogDir(1).folder, catalogDir(1).name) ;
opts = detectImportOptions(catalog_fn);
catalog_data = readtable(catalog_fn, opts) ;

% find row of catalog data that corresponds to current experiment folder
idx = (catalog_data.ExprNum == ExprNum) ;

if sum(idx) ~= 1
    fprintf('Error: could not find unique entry for experiment number \n')
    return
end
% ---------------------------------------------------
%% read info from catalog_data to exprInfoStruct
exprInfoStruct.ExprNum = ExprNum ;

% ------------------
% genotype info
exprInfoStruct.driver = catalog_data.Driver{idx} ;
% effector is a little annoying due to lack of catalog consistency
if ismember('Effector', catalog_data.Properties.VariableNames)
    exprInfoStruct.effector = catalog_data.Effector{idx} ;
elseif contains(exprPath, 'Chrimson')
    exprInfoStruct.effector = 'Chrimson' ;
elseif contains(exprPath, 'Opto Silencing')
    exprInfoStruct.effector = 'GtACR1' ;
else
    exprInfoStruct.effector = [] ;
end

% ------------------------
% pulse(s) info
exprInfoStruct.magPulseStart = catalog_data.magPulseStart(idx) ;
exprInfoStruct.magPulseEnd = catalog_data.magPulseEnd(idx) ;
exprInfoStruct.optoPulseStart = catalog_data.optoPulseStart(idx) ;
exprInfoStruct.optoPulseEnd = catalog_data.optoPulseEnd(idx) ;

% ... make sure that, if we have empty entries, they get stored as nan
pulse_fields = {'magPulseStart', 'magPulseEnd', 'optoPulseStart', ...
    'optoPulseEnd'} ; 
for k = 1:length(pulse_fields)
   if iscell(exprInfoStruct.(pulse_fields{k}))
       exprInfoStruct.(pulse_fields{k}) = nan ; 
   end
end

% -----------------------------
% cine info 
exprInfoStruct.cineFrameStart = catalog_data.cineFrameStart(idx) ; 
exprInfoStruct.cineFrameEnd = catalog_data.cineFrameEnd(idx) ; 

% ---------------------------------
%% save exprInfoStruct?
% it takes a few seconds to calculate, so why do it more than necessary?
if saveFlag
   save(savePathFull, 'exprInfoStruct') 
end

end