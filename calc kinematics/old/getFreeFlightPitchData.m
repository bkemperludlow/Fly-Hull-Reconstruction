rootPath = 'G:\Janelia Flies\kir2.1 flies\Analysis\No Perturbation' ;
savePath = 'G:\Janelia Flies\kir2.1 flies\Analysis\No Perturbation\Compiled Data' ;

saveFlag = false ; 

freeFlightPitchStruct = struct('ExprNum',[],'MovNum',[],'time',[],'bodyPitchAngle',[]) ;

cd(rootPath) ;
movieDir = dir('Expr_8*') ;
Nmovies = length(movieDir) ;

for i = 1:Nmovies 
   folderName = movieDir(i).name ; 
   load([rootPath '\' folderName '\' folderName '_results.mat'])
   set(0, 'ShowHiddenHandles', 'on')
   try
       close 'easyWand 5'
   catch
       disp('No easy wand window')
   end
   ExprNum = str2double(folderName(end-8)) ;
   MovNum = str2double(folderName(end-2:end)) ;
   time = (tin:tout)*(1/8000) ; % in seconds
   axisHats = data.AHat ;
   bodyPitchAngle = (180/pi)*asin(axisHats(:,3)) ; % degrees
   
   freeFlightPitchStruct(i).ExprNum = ExprNum ;
   freeFlightPitchStruct(i).MovNum = MovNum ;
   freeFlightPitchStruct(i).time = time(19:end-18) ; %fucking delta
   freeFlightPitchStruct(i).bodyPitchAngle = bodyPitchAngle ;
   
   if (i == 17) || (i ==31)
       keyboard ;
   end
   
end

if saveFlag
    cd(savePath)
    save freeFlightPitchStruct_Expr_8 freeFlightPitchStruct
end