ExprNum = 1 ;
MovNum = 65 ;
flyType = 1 ; %1 = experimental , 2 = control ;
debugFlag = true ;
plotFlag = true ;

[controller_fit_struct, h_controller, h_summary] = fitRollController(data, ExprNum, MovNum, ...
    flyType, debugFlag, plotFlag) ;

filename = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum)] ;

if (0)
   save([filename  '_controllerFit'], 'controller_fit_struct') ; 
   savefig(h_controller,[filename  '_controllerFit.fig']) ;
   savefig(h_summary, [filename  '_summary.fig']) ;
end
