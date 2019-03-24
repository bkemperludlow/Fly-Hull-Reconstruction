function combineAllMoviesInFolder(folderName, movieNumbers, fileNameString)
% combineAllMoviesInFolder('G:\Janelia Flies\kir2.1 flies\15_08052016', 0:18,'EmptyVector_kir2.1_Expr_15') ;
%combineAllMoviesInFolder('D:\Janelia Flies\tnt round 2\20_01112016\', 47:66,'Expr_20') ;
%combineAllMoviesInFolder('D:\Box Sync Old\VNC Motor Lines\48_12122018\', 70:71,'Expr_48') ;
%prevpath = pwd 

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

for mov = movieNumbers
    try
        if (mov<10)
            zstr = '00' ;
        elseif (mov<100)
            zstr = '0' ;
        else
            zstr='';
        end
        
        disp(['processing movie ' num2str(mov) ]) ;
        xyfile = ['xy_' zstr num2str(mov) '.cin'] ;
        xzfile = ['xz_' zstr num2str(mov) '.cin'] ;
        yzfile = ['yz_' zstr num2str(mov) '.cin'] ;
        
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
        
        outputFileName = [fileNameString '_movie_' zstr num2str(mov)] ;
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