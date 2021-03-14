% -------------------------------------------------------------------------
% function to load data for post processing. want a function to make the
% code for post-processing a bit cleaner.
%
% this will be useful because we often have multiple versions of analysis
% files for the same movie
% -------------------------------------------------------------------------
function [data, analysisOutput, data_out_filename, dataPath, errorFlag] = ...
    hierarchicalLoadData(pathStruct, MovNum, suffixStr_out, ...
     suffixStr_in, loadResultsFlag)
% -------------------------------------------------------------------------
%% clarify inputs
if ~exist('suffixStr_out','var') || isempty(suffixStr_out)
   suffixStr_out = '_cleaned' ;  
end
if ~exist('suffixStr_in','var') || isempty(suffixStr_in)
   suffixStr_in = '_Data_manually_corrected' ;  
end
if ~exist('loadResultsFlag','var') || isempty(loadResultsFlag)
   loadResultsFlag = false ;  
end

% make sure that the suffix strings contain an underscore to keep with
% formating 
if ~strcmp(suffixStr_out(1),'_')
    suffixStr_out = ['_' suffixStr_out] ; 
end
if ~strcmp(suffixStr_in(1),'_')
    suffixStr_in = ['_' suffixStr_in] ; 
end

% boolean value to indicate whether we run into trouble loading data (will
% switch to true if problems arise)
errorFlag = false ;
% ------------------------------------------------------------------------
%% get path to data 
% read experiment number
ExprNum = pathStruct.ExprNum ; 

% find path to folder containing all data files
analysisPath = findMovAnalysisPath(pathStruct, MovNum) ;
if isempty(analysisPath)
    fprintf('Could not find path for expr %d movie %d \n', ExprNum, MovNum)
    data = [] ;
    analysisOutput = [] ;
    data_out_filename = [] ;
    dataPath = [] ; 
    errorFlag = true ;
    return
end

% get prefixes for data (of the form Expr_X_mov_YYY for most, and
% ExprXmovYYY for manually corrected data)
dataPrefix = sprintf('Expr_%d_mov_%03d',ExprNum, MovNum) ; 
dataPrefix_manCorr = sprintf('Expr%dmov%03d',ExprNum, MovNum) ;   

% -----------------------------------------------------------
%% load data (test for multiple different file versions)
% first choice decided by the suffix string entered as funcion input
if contains(suffixStr_in, 'manually_corrected','IgnoreCase',1)
    dataPath = fullfile(analysisPath, ...
        [dataPrefix_manCorr '_Data_manually_corrected.mat']) ; 
else
    dataPath = fullfile(analysisPath, [dataPrefix suffixStr_in '.mat']) ; 
end

% check to see if first choice exists -- if not, come up with ranked list
% of data filenames to try. Note: we may end up trying the same thing
% twice, but not a huge deal
if exist(dataPath,'file') && contains(suffixStr_in, '_results')
    analysisOutput = load(dataPath) ;
    data = analysisOutput.data ;
elseif exist(dataPath,'file')
    data = importdata(dataPath) ;
else
    % first choice: manuallyCorrected
    dataPath1 = fullfile(analysisPath, ...
        [dataPrefix_manCorr '_Data_manually_corrected.mat']) ;
    % second choice: cleaned
    dataPath2 = fullfile(analysisPath, [dataPrefix '_cleaned.mat']) ;
    % third choice: test
    dataPath3 = fullfile(analysisPath, [dataPrefix '_test.mat']) ;
    % third choice: results
    dataPath4 = fullfile(analysisPath, [dataPrefix '_results.mat']) ;
    
    % check the existence of these file paths. if none exist, we'll have to
    % use full analysis output
    if exist(dataPath1,'file')
        data = importdata(dataPath1) ;
        dataPath = dataPath1 ; 
    elseif exist(dataPath2,'file')
        data = importdata(dataPath2) ;
        dataPath = dataPath2 ;
    elseif exist(dataPath3,'file')
        data = importdata(dataPath3) ;
        dataPath = dataPath3 ;
    elseif exist(dataPath4,'file')
        analysisOutput = load(dataPath4) ;
        data = analysisOutput.data ;
        dataPath = dataPath4 ;
    else
        fprintf('Could not find data file for expr %d movie %d \n', ...
            ExprNum, MovNum)
        data = [] ;
        analysisOutput = [] ;
        data_out_filename = [] ;
        dataPath = [] ; 
        errorFlag = true ;
        return
    end
end

% -----------------------------------------------------------
%% load "analysisOutput" (?)
if loadResultsFlag
   resultsPath = fullfile(analysisPath, [dataPrefix '_results.mat']) ;
   if exist(resultsPath, 'file')
    analysisOutput = load(resultsPath) ; 
   else
       analysisOutput = [] ; 
       errorFlag = true ; 
   end
else
   analysisOutput = [] ; 
end

% -----------------------------------------------------------
%% generate output filename
if contains(suffixStr_out, 'manually_corrected','IgnoreCase',1)
    data_out_filename = fullfile(analysisPath, [dataPrefix_manCorr ...
        suffixStr_out '.mat']) ;
else
    data_out_filename = fullfile(analysisPath, [dataPrefix ...
        suffixStr_out '.mat']) ;
end


end