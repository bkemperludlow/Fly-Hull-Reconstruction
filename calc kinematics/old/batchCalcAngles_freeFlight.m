rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Right' ;


cd(rootPath) ;
movieDir = dir('Expr*') ;
Nmovies = length(movieDir) ;

for i = 1:length(movieDir)
   folderName = movieDir(i).name ; 
   datapath = [rootPath '\' folderName '\'] ;
   exprNum = str2double(folderName(end-9:end-8)) ;
   if isnan(exprNum)
       exprNum = str2double(folderName(end-8)) ;
   end
   movNum = str2double(folderName(end-2:end)) ;
   
   datafilename = [datapath '\' folderName '_test.mat'] ;
    
   cd(datapath)
   
   if exist(datafilename,'file') == 2
       continue ;
   else 
       try
           load([rootPath '\' folderName '\' folderName '_results.mat'])
       catch
           continue ; 
       end
       set(0, 'ShowHiddenHandles', 'on')
       try
           close 'easyWand 5'
       catch
           disp('No easy wand window')
       end
       
       [rhoTimes, rollVectors] = estimateRollVector(data) ;
       data.rhoTimes = rhoTimes ;
       data.rollVectors = rollVectors ;
       
       run quick_and_dirty
       
       save(datafilename, 'data')
       %keyboard ;
       close all
   end
end