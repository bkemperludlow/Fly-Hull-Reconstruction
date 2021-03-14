function combineAllMoviesInFolder(folderName, movieNumbers, ...
    fileNameString, movFileExt)
% combineAllMoviesInFolder('D:\Box Sync Old\Opto Silencing\43_24022020', 6:44,'Expr43') ;
%combineAllMoviesInFolder('D:\Janelia Flies\tnt round 2\20_01112016\', 47:66,'Expr_20') ;
%combineAllMoviesInFolder('D:\Box Sync Old\Opto Silencing\38_28012020\', 44:62,'Expr_38') ;
%{
combineAllMoviesInFolder('D:\Box Sync Old\Opto Silencing\49_03112020\');
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
    movFileExt = '.cin' ; 
end

% ----------------------------------------------------------------------    
cd(folderName) ;

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
        
        outputFileName = [fileNameString '_movie_' num2str(mov,'%03d')] ;
        fps = 30 ;
        realFrameRate = 8000 ;
   
        
        combineMoviesFromCin(metaData, outputFileName,fps, realFrameRate) ;
        
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