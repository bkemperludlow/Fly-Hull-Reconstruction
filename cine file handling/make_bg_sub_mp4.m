% -------------------------------------------------------------------------
% function to generate mp4 with background removed and colored square to
% indicate perturbation (opto or stim)
% -------------------------------------------------------------------------
function [] = make_bg_sub_mp4(rootPath, MovNum, pulseDuration, ...
    patchColor, output_fps)
% -----------------------------
%% inputs and params
if ~exist('pulseDuration','var') || isempty(pulseDuration)
    pulseDuration = 50 ; % in milliseconds
end
if ~exist('patchColor','var') || isempty(patchColor)
    patchColor = [] ; % color for 'insertShape'. If empty, won't draw patch
    % NB: can be 'blue' | 'green' | 'red' | 'cyan' | 'magenta' | 'black' | 'black' | 'white'
end
if ~exist('output_fps','var') || isempty(output_fps)
    output_fps = 50 ; % frames per second for output mp4
end
% -------------------------------------------------------------------------
%% get path for movie 1) cine files and 2) analysis
pathStruct = generatePathStruct(rootPath) ; 
movAnalysisPath = findMovAnalysisPath(pathStruct, MovNum) ; 
ExprNum = pathStruct.ExprNum ;

% directory containing data info
dataDir = dir(fullfile(movAnalysisPath, '*_results.mat')) ;

% move to root directory so we don't write random image files in the code
% base
prevpath = pwd ;
cd(rootPath) ;

% prefix name for movies
fileNameString = ['Expr_' num2str(ExprNum)] ;
outputFileName = [fileNameString '_movie_' num2str(MovNum,'%03d') '_bgSub'] ;

% see if analysis data exist
if length(dataDir) ~= 1
    fprintf('Could not locate analysis  for Expr %d, Mov %d--quitting \n',...
        ExprNum, MovNum)
    return
end
data_in = importdata(fullfile(dataDir.folder,dataDir.name)) ;
allBG = data_in.allBG ;
realFrameRate = data_in.data.params.fps ;
%--------------------------------------------------------------------------
%% initialize cine data 
s.height = size(allBG,2) ;
s.width  = size(allBG,3) ;
s.exists = 1 ;
s.firstImage = -1 ;
s.lastImage = -1 ;
s.filename = ' ' ;
s.cindata = [] ;

metaData = repmat(s,1,3) ;

% camera indices
XZ = 1 ;
XY = 2 ;
YZ = 3 ;

%--------------------------------------------------------------------------
%% load cine files
xyfile = ['xy_' num2str(MovNum,'%03d') '.cin'] ;
xzfile = ['xz_' num2str(MovNum,'%03d') '.cin'] ;
yzfile = ['yz_' num2str(MovNum,'%03d') '.cin'] ;

try
    xydat = getCinMetaData(xyfile) ;
    s.firstImage = xydat.firstImage ;
    s.lastImage  = xydat.lastImage ;
    s.filename   = xyfile ;
    s.exists     = 1 ;
    s.cindata = myOpenCinFile(xyfile) ;
catch
    s.exists = 0 ;
    s.cindata = [] ;
end
metaData(XY) = s ; %#ok<AGROW>

try
    xzdat = getCinMetaData(xzfile) ;
    s.firstImage = xzdat.firstImage ;
    s.lastImage  = xzdat.lastImage ;
    s.filename   = xzfile ;
    s.exists     = 1 ;
    s.cindata = myOpenCinFile(xzfile) ;
catch
    s.exists = 0 ;
    s.cindata = [] ;
end
metaData(XZ) = s ; %#ok<AGROW>

try
    yzdat = getCinMetaData(yzfile) ;
    s.firstImage = yzdat.firstImage ;
    s.lastImage  = yzdat.lastImage ;
    s.filename   = yzfile ;
    s.cindata = myOpenCinFile(yzfile) ;
    s.exists = 1 ;
catch
    s.exists = 0 ;
    s.cindata = [] ;
end
metaData(YZ) = s ; %#ok<AGROW>

% ---------------------------------------------------------------------
%% run function to combine movies into mp4 in chunks
combineMoviesFromCin_bgSub(metaData, outputFileName,output_fps, ...
    realFrameRate, allBG, pulseDuration, patchColor ) ;

% -------------------------------------------------------
%% close cine files
if (metaData(XY).exists)
    myCloseCinFile(metaData(XY).cindata) ;
end

if (metaData(XZ).exists)
    myCloseCinFile(metaData(XZ).cindata) ;
end

if (metaData(YZ).exists)
    myCloseCinFile(metaData(YZ).cindata) ;
end

% ------------------------------------------------------
%% return to original directory location
cd(prevpath)
end