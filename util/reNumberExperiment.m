% assumes data has already been moved to correct folder (i.e. rootPath is
% the folder that you want the data in)
rootPath = 'D:\Box Sync\VNC Sensory Lines\11_01082018\' ;
ExprNumOld = 34 ;

pathSplit = strsplit(rootPath,'\') ;
folderSplit = strsplit(pathSplit{end-1},'_') ;
ExprNumNew = str2double(folderSplit{1}) ;

prevPath = pwd ;

%% get folders with incorrect numbering
analysisPath = [rootPath 'Analysis\'] ;
analysisDir = dir([rootPath 'Expr_' num2str(ExprNumOld) '*']) ;
unsortedPath = [analysisPath 'Unsorted\'] ;
unsortedDir = dir([unsortedPath 'Expr_' num2str(ExprNumOld) '*']) ;

%% rename folders
cd(rootPath)
for i = 1:length(analysisDir)
    folderName = analysisDir(i).name ;
    folderName_new = ['Expr_' num2str(ExprNumNew) '_mov_' folderName(end-2:end)] ;
    
    cmd_str=['ren ' folderName ' ' folderName_new];
    
    dos (cmd_str);
    
end

cd(unsortedPath)
for i = 1:length(unsortedDir)
    folderName = unsortedDir(i).name ;
    folderName_new = ['Expr_' num2str(ExprNumNew) '_mov_' folderName(end-2:end)] ;
    
    cmd_str=['ren ' folderName ' ' folderName_new];
    
    dos (cmd_str);
end

%% rename data files
dataFilesDir = rdir([unsortedPath 'Expr_' num2str(ExprNumNew)...
    '*\Expr_' num2str(ExprNumOld) '*']) ;

for j = 1:length(dataFilesDir)
    [filePath, fileName, fileExt] = fileparts(dataFilesDir(j).name) ;
    fileName_split = strsplit(fileName,'_') ;
    
    fileName_new = [fileName_split{1} '_' num2str(ExprNumNew)] ;
    for k = 3:length(fileName_split)
        fileName_new = [fileName_new '_' fileName_split{k} ] ;
    end
    
    cd(filePath)
    cmd_str = ['ren ' fileName fileExt ' ' fileName_new fileExt] ;
    dos (cmd_str);
end

cd(prevPath)
