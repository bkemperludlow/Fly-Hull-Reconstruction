function fname = getDataFileName (exprNumber, movieNumber, PC)
% returns the file name of the "data" structure corresponding to
% the specified movie and experiment
% PC is
%   1 = lab computer (default)
%   2 = laptop
%
% this function also sets the path.

fname = -1 ; % value returned upon error

if (~exist('PC','var'))
    PC = 1;
end

switch PC
    case 1
        addpath F:\Roll\MatlabCode\TRACKER\analysisCode\
        addpath F:\Roll\MatlabCode\TRACKER\
        addpath F:\Roll\MatlabCode\TRACKER\fdaM
        basePath = ['F:\Roll\MatlabCode\TRACKER\Expr' num2str(exprNumber) '\'];
    case 2
        addpath C:\Tsevi\Roll\MatlabCode\TRACKER\analysisCode\
        addpath C:\Tsevi\Roll\MatlabCode\TRACKER\
        addpath C:\Tsevi\Roll\MatlabCode\TRACKER\fdaM
        basePath = ['C:\Tsevi\Roll\MatlabCode\TRACKER\Expr' num2str(exprNumber) '\'] ;
    otherwise
        disp ('Invalid value for PC. Aborting.') ;
        return
end

folderName = ['expr' num2str(exprNumber) 'mov' num2str(movieNumber) '\'];

combinedPath = [basePath folderName ] ;

switch (exprNumber)
    case 22 % experiment no. 22
        switch(movieNumber)
            case 34
                outputFilename = 'mov34data_Oct2011' ;
            case 26
                outputFilename = 'mov26data_Oct2011' ;
            case 4
                outputFilename = 'mov4data_revisited_April_2014' ; % 'mov4data' ;
            case 23
                %outputFilename = 'mov23data_Oct2011b' ;
                outputFilename = 'Oct2012\mov23data_revisited_Oct_2012' ;
            case 5
                outputFilename = 'mov5_user_corrected' ;
            case 32
                outputFilename = 'mov32_user_corrected' ;
            otherwise
                disp('invalid movie number. aborting.')
                return
        end
        
    case 27 % experiment no. 27
        switch(movieNumber)
            case 22
                outputFilename = 'mov22_user_corrected' ;
            case 58
                outputFilename = 'mov58_user_corrected' ;
            case 102
                outputFilename = 'mov102_user_corrected' ;
            case 117
                outputFilename = 'mov117_user_corrected' ;
            case 48
                outputFilename = 'mov48data' ;
            case 32 
                outputFilename = 'expr27mov32data' ;
            case 46
                outputFilename = 'expr27mov46data' ;
            otherwise
                disp('invalid movie number. aborting.')
                return
        end
    case 28
        switch(movieNumber)
            case 47
                outputFilename = 'expr28mov47dataCombined' ; % 'expr28mov47data';
            case 56
                outputFilename = 'expr28mov56data' ;
            case 63
                outputFilename = 'expr28mov63data' ;
            case 71
                outputFilename = 'expr28mov71data_trimmed' ; % 'expr28mov71data' ;
            case 73
                outputFilename = 'expr28mov73data_trimmed' ;
            case 168
                outputFilename = 'mov168_user_corrected' ;
            case 190
                outputFilename = 'mov190_user_corrected_trimmed' ;
            otherwise
                disp('invalid movie number. aborting.')
                return
        end
    otherwise
        disp('invalid experiment number. aborting.') ;
        return
        
end
fname = [ combinedPath outputFilename '.mat'] ;

end
