function combineAllMoviesInFolder(folderName, movieNumbers, ...
    fileNameString, movFileExt, realFrameRate)

%combineAllMoviesInFolder('D:\Box Sync Old\Opto Silencing\38_28012020\', 44:62,'Expr_38') ;
%{
combineAllMoviesInFolder('D:\Box Sync Old\Mechanical Perturbation\04_30072024\');

combineAllMoviesInFolder('D:\Box Sync Old\Opto Mechanical\195_27072024\')
combineAllMoviesInFolder('D:\Box Sync Old\Haltere Experiment\91_02082024\')
combineAllMoviesInFolder('Y:\Optomechanical\155_09112023\')


Y:\Optomechanical\150_31102023
Y:\Optomechanical\152_04112023
Y:\Optomechanical\151_02112023
Y:\Optomechanical\153_05112023

combineAllMoviesInFolder('D:\Box Sync Old\Antenna test\01_06082024')

%}
%prevpath = pwd 
% -------------------------------------------------------------------------
%% inputs and params
% movie numbers
if ~exist('movieNumbers','var') || isempty(movieNumbers)
    % if we don't get input, just assume we should do mp4s for all movies
    xydir = dir(fullfile(folderName, 'xy_*.cin*')) ; 
    movNumExpression = '(?<camName>[xyz]{2})_(?<movieNum>\d+)' ; 
    movNumStrs = arrayfun(@(x) regexp(x.name, movNumExpression,'names'), ...
        xydir, 'UniformOutput',0) ; 
    movieNumbers = cellfun(@(y) str2double(y.movieNum), movNumStrs) ; 
    
    % since we're looping over "movieNumbers" need to make sure it's a row
    % vector
    if size(movieNumbers,1) > size(movieNumbers,2)
       movieNumbers = movieNumbers' ;  
    end
end
% input for fileNameString
if ~exist('fileNameString','var') || isempty(fileNameString)
    % pattern to search for in Experiment folder string
    exprFolderExpression = '(?<exprNum>\d+)_(?<datenum>\d+)' ;
    
    % isolate experiment folder (XX_DDMMYYYY)
    folderNameSplit = strsplit(folderName,'\') ; 
    empty_idx = cellfun(@(y) isempty(y), folderNameSplit) ; 
    folderNameSplit = folderNameSplit(~empty_idx) ; 
    exprFolder = folderNameSplit{end} ; 
    
    % find experiment number from folder name
    reg_exp_out = regexp(exprFolder, exprFolderExpression ,'names') ; 
    ExprNum = reg_exp_out.exprNum ; 
    
    % create default fileNameString from experiment number
    fileNameString = strjoin({'Expr', num2str(ExprNum)}, '_') ; 
end
if ~exist('movFileExt','var') || isempty(movFileExt)
    % get file type (.cin vs .cine) if we don't have it
    xydir = dir(fullfile(folderName, 'xy_*.cin*')) ; 
    extSearchStr = '(?<fileName>\w+).(?<fileExt>\w+)' ; 
    fnPartsList = arrayfun(@(x) regexp(x.name,extSearchStr,'names'), ...
        xydir, 'UniformOutput',false) ; 
    fileExtList = cellfun(@(y) y.fileExt, fnPartsList, 'UniformOutput',0) ; 
    fileExtListUnique = unique(fileExtList) ;
    
    % make sure we only have one file type in folder
    if length(fileExtListUnique) ~= 1
       fprintf('Warning: multiple phantom file types in %s \n', folderName)
       keyboard
    end
    movFileExt = strcat('.',fileExtListUnique{1}) ; 
end
if ~exist('realFrameRate','var') || isempty(realFrameRate)
    % frame capture rate of cameras
    realFrameRate = 8000 ;
end

% ----------------------------------------------------------------------    
cd(folderName) ;

% output frame rate (for mp4)
fps = 30 ;

% default info for cines (but can be overwritten)
s.height = 512 ;
s.width  = 512 ; 
s.exists = 1 ;
s.firstImage = -1 ;
s.lastImage = -1 ;
s.filename = ' ' ;
s.cindata = [] ;

metaData = repmat(s,1,3) ;

XZ = 1 ;
XY = 2 ;
YZ = 3 ;

% -------------------------------------------------------------------------
%% loop over movies
for mov = movieNumbers
    try
        disp(['processing movie ' num2str(mov) ]) ;
        xyfile = ['xy_' num2str(mov,'%03d') movFileExt] ;
        xzfile = ['xz_' num2str(mov,'%03d') movFileExt] ;
        yzfile = ['yz_' num2str(mov,'%03d') movFileExt] ;
        
        try
            xydat = getCinMetaData(xyfile) ;
            s.firstImage = xydat.firstImage ;
            s.lastImage  = xydat.lastImage ;
            s.height     = xydat.height ; 
            s.width      = xydat.width ; 
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
            s.height     = xzdat.height ; 
            s.width      = xzdat.width ; 
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
            s.height     = yzdat.height ; 
            s.width      = yzdat.width ; 
            s.filename   = yzfile ;
            s.cindata = myOpenCinFile(yzfile) ;
            s.exists = 1 ;
        catch
            s.exists = 0 ;
            s.cindata = [] ;
        end   
        metaData(YZ) = s ; %#ok<AGROW>
        
        % initialize output file
        outputFileName = [fileNameString '_movie_' num2str(mov,'%03d')] ;
        
        % combine videos
%         combineMoviesFromCin(metaData, outputFileName,fps, realFrameRate) ;
        combineMoviesFromCinSam(metaData, outputFileName,fps, realFrameRate) ;
        
        % close cine objects
        if (metaData(XY).exists)
            myCloseCinFile(metaData(XY).cindata) ;
        end
        
        if (metaData(XZ).exists)
            myCloseCinFile(metaData(XZ).cindata) ;
        end
        
        if (metaData(YZ).exists)
            myCloseCinFile(metaData(YZ).cindata) ;
        end
        
        
    catch ERR
        disp(ERR) ;
        disp('cannot process movie...')
        keyboard ;
    end
end
%cd(prevpath) ;

%=========================================
% move mp4 files to folder post conversion
mp4Path = [folderName '\mp4\'] ; 
if exist(mp4Path,'dir') > 0 
    movefile([folderName '\*.mp4'], mp4Path)
end

end