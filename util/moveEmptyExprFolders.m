%--------------------------------------------------------------------------
% when analysis fails to run, folders pile up in root analysis folder. this
% moves them so i don't have to wonder which folders have data in them
%--------------------------------------------------------------------------
function [] = moveEmptyExprFolders(rootPath)
%{
rootPath = 'D:\Box Sync\VNC MN Chrimson\12_23072018' ; 
moveEmptyExprFolders(rootPath)
%}
analysisPath = fullfile(rootPath, 'Analysis') ; 
failPath = fullfile(analysisPath, 'Failed Analysis') ; 
unsortPath = fullfile(analysisPath, 'Unsorted') ; 

MIN_MEMORY = 1024 ; % bytes

mkdir(failPath) ; 
rootPathDir = rdir([analysisPath '\Expr*']) ;  

for i = 1:length(rootPathDir)
   dataFolder = rootPathDir(i).name ;
   dataFolderDir = dir(dataFolder) ; 
   
   if (sum([dataFolderDir(:).bytes]) < MIN_MEMORY)
       movefile(dataFolder, failPath) ; 
   else
       try 
           movefile(dataFolder, unsortPath) ; 
       catch
           continue
       end
   end
    
end
end