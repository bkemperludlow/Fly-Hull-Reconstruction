folderName = 'G:\Janelia Flies\kir2.1 flies round 2\08_07072016\' ;
dataPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\No Perturbation\Expr_8_mov_023\' ;

mov = 23 ;
patch_color = 'red' ; % or e.g. 'green' for uas-gtACR1
prevpath = pwd ;
cd(folderName) ;

optoFlag = false ;

fileNameString = 'Expr_8' ;
dataDir = dir([dataPath '*_results.mat']) ;

outputFileName = [fileNameString '_movie_' num2str(mov,'%03d')] ;
fps = 50 ;
realFrameRate = 8000 ;


if length(dataDir) ~= 1
    keyboard
end
data_in = importdata(fullfile(dataDir.folder,dataDir.name)) ;
allBG = data_in.allBG ;
%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------

xyfile = ['xy_' num2str(mov,'%03d') '.cin'] ;
xzfile = ['xz_' num2str(mov,'%03d') '.cin'] ;
yzfile = ['yz_' num2str(mov,'%03d') '.cin'] ;

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

if optoFlag
    combineMoviesFromCin_bgSub_v2(metaData, outputFileName,fps, realFrameRate,...
        allBG,patch_color) ;
else
    combineMoviesFromCin_bgSub(metaData, outputFileName,fps, realFrameRate,...
        allBG) ;
end

if (metaData(XY).exists)
    myCloseCinFile(metaData(XY).cindata) ;
end

if (metaData(XZ).exists)
    myCloseCinFile(metaData(XZ).cindata) ;
end

if (metaData(YZ).exists)
    myCloseCinFile(metaData(YZ).cindata) ;
end

cd(prevpath)
