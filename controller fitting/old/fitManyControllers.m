%Run controllerFit_v2 many times

Expr = [1];
Mov = [31] ;
PitchType = 'up' ;
saveFlag = false ;

if strcmp(PitchType,'up')
    savePathRoot = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Up\' ;
elseif strcmp(PitchType,'down')
    savePathRoot = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Down\' ;
end

runSize = length(Expr) ; 

resultsArray = zeros(runSize, 8) ; 
resultsArray(:,1) = Expr' ;
resultsArray(:,2) = Mov' ;

CIArray = zeros(runSize,3) ;  
%{
resultsArray_NoK = zeros(runSize, 9) ; 
resultsArray_NoK(:,1) = Expr' ;
resultsArray_NoK(:,2) = Mov' ;

resultsArray_NoKi = zeros(runSize, 9) ; 
resultsArray_NoKi(:,1) = Expr' ;
resultsArray_NoKi(:,2) = Mov' ;

resultsArray_NoKp = zeros(runSize, 9) ; 
resultsArray_NoKp(:,1) = Expr' ;
resultsArray_NoKp(:,2) = Mov' ;

resultsArray_NoKp_NoK = zeros(runSize, 9) ; 
resultsArray_NoKp_NoK(:,1) = Expr' ;
resultsArray_NoKp_NoK(:,2) = Mov' ;

resultsArray_SameT = zeros(runSize, 9) ; 
resultsArray_SameT(:,1) = Expr' ;
resultsArray_SameT(:,2) = Mov' ;

resultsArray_SameT_NoK = zeros(runSize, 9) ; 
resultsArray_SameT_NoK(:,1) = Expr' ;
resultsArray_SameT_NoK(:,2) = Mov' ;
%}
for i = 1:runSize
    %{
    [results_full, results_NoK, results_NoKi, results_SameT,...
    results_SameT_NoK, results_NoKp, results_NoKp_NoK, PitchEstErr] ...
        = controllerFit_v2(Expr(i),Mov(i),PitchType,true,false,false,true) ;
    %}
    [results_full, CI, PitchEstErr] = controllerFit_v2(Expr(i),Mov(i),PitchType,false,...
        true,true,true) ;
    
    resultsArray(i,3:7) = results_full ;
    resultsArray(i,8) = PitchEstErr ;
    
    CIArray(i,:) = CI ; 
    
    %{
    resultsArray_NoK(i,3:8) = results_NoK ;
    resultsArray_NoK(i,9) = PitchEstErr ;
    
    resultsArray_NoKi(i,3:8) = results_NoKi ;
    resultsArray_NoKi(i,9) = PitchEstErr ;
    
    resultsArray_NoKp(i,3:8) = results_NoKp ;
    resultsArray_NoKp(i,9) = PitchEstErr ;
    
    resultsArray_NoKp_NoK(i,3:8) = results_NoKp_NoK ;
    resultsArray_NoKp_NoK(i,9) = PitchEstErr ;
    
    resultsArray_SameT(i,3:8) = results_SameT ;
    resultsArray_SameT(i,9) = PitchEstErr ;
    
    resultsArray_SameT_NoK(i,3:8) = results_SameT_NoK ;
    resultsArray_SameT_NoK(i,9) = PitchEstErr ;
    %}
    %{
    if saveFlag
        if Mov(i) < 10 
            zstr = '00' ;
        elseif Mov(i) < 100
            zstr = '0' ;
        else 
            zstr = '' ;
        end
        
        savePath = strcat(savePathRoot,'Expr_',num2str(Expr(i)),...
            '_mov_',zstr,num2str(Mov(i)),'\control fit') ;
        try 
            cd(savePath)
        catch
            mkdir('control fit')
            cd(savePath)
        end
        save controlFits resultsArray resultsArray_NoK resultsArray_NoKi resultsArray_NoKp...
            resultsArray_NoKp_NoK resultsArray_SameT resultsArray_SameT_NoK
        
        fig1name = ['Expr ' num2str(Expr(i)) ' Mov ' num2str(Mov(i)) ' Contributions'] ;
        fh1 = findobj( 'Type', 'Figure', 'Name', fig1name );
        fig2name = ['Expr ' num2str(Expr(i)) ' Mov ' num2str(Mov(i)) ' Contributions 2'] ;
        fh2 = findobj( 'Type', 'Figure', 'Name', fig2name );
        
        figure(fh1(end)) ;
        savefig 'controlFit1'
        figure(fh2(end)) ;
        savefig 'controlFit2'
    end
    %}
end
%cd('F:\luca\Analysis\metadata\Controller\24-1-2015')
%save controlFits resultsArray CIArray
%{
save controlFits resultsArray resultsArray_NoK resultsArray_NoKi resultsArray_NoKp...
    resultsArray_NoKp_NoK resultsArray_SameT resultsArray_SameT_NoK
%}