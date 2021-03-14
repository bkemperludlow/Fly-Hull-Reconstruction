function [] = cleanUpAnalysesFun(pathToWatch, analysisType, MovNum, ...
    clustFlag, largePertFlag, removeLegsFlag, alignBBoxFlag)
% -------------------------------------------------------------------------
% FUNCTION to perform different types of post-processing to fly data
%
%  analysisType options: 'new_recon' | 'full' | 'full_new' | 'just_angles'
%   | 'stop_wings' | 'clean_wings'  | 'refine_wing_vecs' | 
%    'redo_hull_analysis' | 'swap_wings' | 'extreme_roll' | 'refine_chord'
% -------------------------------------------------------------------------
%% input parameters/paths
% -----------------------------
% set path for experiment folder wherein analysis is still needed
if ~exist('pathToWatch','var') || isempty(pathToWatch)
    pathToWatch = 'D:\Box Sync Old\Opto Silencing\46_23102020\' ; 
end
% do full analysis, new reconstruction method, just angles, or other?
if ~exist('analysisType','var') || isempty(analysisType)
    analysisType = 'clean_wings' ; % 'extreme_roll' ; %'clean_wings' ;
end
% range of movies to analyze in given experiment?
if ~exist('MovNum','var') || isempty(MovNum)
    MovNum = (122:134) ; % 
end
if ~exist('clustFlag','var') || isempty(clustFlag)
    clustFlag = false ; % false % which version of analysis script to run
end
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false  ; % is it a large perturbation?
end
if ~exist('removeLegsFlag','var') || isempty(removeLegsFlag)
    removeLegsFlag = true ; % try to remove legs in binary threshold?
end
if ~exist('alignBBoxFlag','var') || isempty(alignBBoxFlag)
    alignBBoxFlag = false ; % try to align images to avoid clipping?
end

% get experiment number from folder name--too lazy to re-enter each time
pathStruct = generatePathStruct(pathToWatch) ;
pathSplit = strsplit(pathToWatch,'\') ;
folderSplit = strsplit(pathSplit{end-1},'_') ;
ExprNum = str2double(folderSplit{1}) ;

% set movie numbers that need to be analyzed
Nmovies = length(MovNum) ;

% run analysis of movies
for k = 1:Nmovies
    %tic
    fprintf('\n Starting analysis on movie number: %d \n \n', MovNum(k))
    switch analysisType
        %% --------------------------------------------------
        % just calculate angles from extant reconstruction
        % --------------------------------------------------
        case 'just_angles'
            analysisPath = findMovAnalysisPath(pathStruct, MovNum(k)) ;
            if isempty(analysisPath)
                fprintf('Could not find path for movie %d \n', MovNum(k))
                continue
            end
            analysisPath_split = strsplit(analysisPath,'\') ;
            savePath = strjoin({analysisPath_split{1:end-1}},'\') ;
            try
                estimateFlyAngles(ExprNum, MovNum(k), savePath,...
                    largePertFlag) ;
            catch
                fprintf('Failed to analyze movie %d \n', MovNum(k))
            end
        %% ----------------------------------------------------------
        % perform full fly analysis from cine files to kinematics
        % ----------------------------------------------------------
        case 'full'
            flyAnalysisMain(MovNum(k), ExprNum, pathStruct, clustFlag, ...
                largePertFlag, removeLegsFlag) ;
        %% ----------------------------------------------------------------
        % perform reconstruction using image shift method to avoid voxel
        % clipping
        % ----------------------------------------------------------------
        case 'new_recon'
            data = redoHullReconstruction(MovNum(k), ExprNum, ...
                pathStruct, alignBBoxFlag) ;
        %% ---------------------------------------------------------------
        % perform reconstruction using modified binaryThreshold--should
        % take care of issues of wing identification when wings stop
        % flapping
        % ---------------------------------------------------------------
        case 'full_new'
            flyAnalysisMain(MovNum(k), ExprNum, pathStruct, clustFlag,...
                largePertFlag, removeLegsFlag, false, true) ;
        %% ---------------------------------------------------------------
        % perform reconstruction using modified binaryThreshold--should
        % take care of issues of wing identification when wings stop
        % flapping
        % ---------------------------------------------------------------
        case 'stop_wings'
            flyAnalysisMain(MovNum(k), ExprNum, pathStruct, clustFlag,...
                largePertFlag, removeLegsFlag, true) ;
        %% ---------------------------------------------------------------
        % use extant data structure and clean up wing voxels using
        % clustering/swapping
        % ---------------------------------------------------------------
        case 'clean_wings'
            % ---------------------------------------------------
            % load data
            suffixStr_out = '_cleaned' ; 
            suffixStr_in = '_new' ; %'_test' ; 
            loadResultsFlag = true ; 
            [data_in, analysisOutput, data_out_filename, dataPath,...
                errorFlag] = hierarchicalLoadData(pathStruct, MovNum(k), ...
                suffixStr_out, suffixStr_in, loadResultsFlag) ; 
            
            if errorFlag
                %fprintf('Could not find path for movie %d \n', MovNum(k))
                continue
            end
            
            % --------------------------------------------
            % read in info from "results" file
            dlt = analysisOutput.dlt_matrix ; 
            all_fly_bw = analysisOutput.all_fly_bw ; 
            body_only_bw = analysisOutput.body_only_bw ; 
            order = [2, 1, 3] ;
            clear analysisOutput
            
            % --------------------------------
            % run function to clean voxels
            data = cleanUpWingVoxels(data_in, all_fly_bw, ...
                body_only_bw, dlt, order) ;
            
            % -------------------------------------------
            % check to see if wings need to be swapped
            fprintf('Checking L<->R wing swaps for movie %d... \n', MovNum(k))
            swapFlag = true(data.Nimages,1) ; 
            cc = 0 ;
            while (sum(swapFlag) > 1) && (cc < 20)
                cc = cc + 1 ; 
                [data, swapFlag] = checkWingSwap(data, largePertFlag) ;
            end
            
            % ---------------------------
            % re-calculate angles
            fprintf('Re-calculating angles for movie %d... \n', MovNum(k))
%             [rhoTimes, rollVectors] = estimateRollVector(data) ;
            [rhoTimes, rollVectors] = estimateRollVector(data,largePertFlag, true) ;
            data.rhoTimes = rhoTimes ;
            data.rollVectors = rollVectors ;
            
            save(data_out_filename, 'data')
            data = calcAnglesMain(data_out_filename, largePertFlag, ...
                true, true, suffixStr_out) ;
            close all
        %% ---------------------------------------------------------------
        % Iterative refinement of wing vectors
        % ---------------------------------------------------------------
        case 'refine_wing_vecs'
            % -------------------------
            % load data 
            suffixStr_out = '_iterRefine' ; 
            suffixStr_in  = 'cleaned' ; 
            [data_in, ~, data_out_filename, ~, errorFlag] = ...
                hierarchicalLoadData(pathStruct, MovNum(k), ...
                 suffixStr_out, suffixStr_in) ; 
            
             if errorFlag
                %fprintf('Could not find path for movie %d \n', MovNum(k))
                continue
             end
            % ----------------------------
            % run function
            data = refineWingVecIter_v2(data_in, [], true) ;
            
            % ----------------------------
            % calculate angles
%             [rhoTimes, rollVectors] = estimateRollVector(data) ;
%             data.rhoTimes = rhoTimes ;
%             data.rollVectors = rollVectors ;
            
            save(data_out_filename, 'data')
            data = calcAnglesMain(data_out_filename, largePertFlag, ...
                true, true, suffixStr_out) ;
            %close all
        
        %% ---------------------------------------------------------------
        % Re-run hullAnalysis_mk*.m 
        % ---------------------------------------------------------------
        case 'redo_hull_analysis'
            % ----------------------------
            % load data
            suffixStr_out = '_test' ;
            suffixStr_in = '_results' ; 
            loadResultsFlag = true ; 
            [~, analysisOutput, data_out_filename, ~, errorFlag] = ...
                hierarchicalLoadData(pathStruct, MovNum(k), ...
                 suffixStr_out, suffixStr_in, loadResultsFlag) ; 
            
            if errorFlag
                %fprintf('Could not find path for movie %d \n', MovNum(k))
                continue
            end
            % ----------------------------
            % re run hull analysis
            data = redoHullAnalysis(analysisOutput) ; 
            
            % ---------------------------
            % re-calculate angles
            [rhoTimes, rollVectors] = estimateRollVector(data,largePertFlag, true) ;
            data.rhoTimes = rhoTimes ;
            data.rollVectors = rollVectors ;
            
            save(data_out_filename, 'data')
            data = calcAnglesMain(data_out_filename, largePertFlag, ...
                true, true) ;
            close all
            
        %% ---------------------------------------------------------------
        % Just see if left and right wings need to be swapped
        % ---------------------------------------------------------------
        case 'swap_wings' 
            % -------------------------
            % load data 
            suffixStr_out = '_cleaned' ;  
            suffixStr_in  = '_cleaned' ; 
            [data, ~, data_out_filename, ~, errorFlag] = ...
                hierarchicalLoadData(pathStruct, MovNum(k), ...
                 suffixStr_out, suffixStr_in) ; 
            
             if errorFlag
                %fprintf('Could not find path for movie %d \n', MovNum(k))
                continue
             end
             
            % -------------------------------------------
            % check to see if wings need to be swapped
            fprintf('Checking L<->R wing swaps for movie %d... \n', MovNum(k))
            swapFlag = true(data.Nimages,1) ; 
            cc = 0 ;
            while (sum(swapFlag) > 1) && (cc < 20)
                cc = cc + 1 ; 
                [data, swapFlag] = checkWingSwap(data, largePertFlag) ;
            end
            
            % ---------------------------
            % re-calculate angles
            fprintf('Re-calculating angles for movie %d... \n', MovNum(k))
%             [rhoTimes, rollVectors] = estimateRollVector(data) ;
            [rhoTimes, rollVectors] = estimateRollVector(data,largePertFlag, true) ;
            data.rhoTimes = rhoTimes ;
            data.rollVectors = rollVectors ;
            
            save(data_out_filename, 'data')
            data = calcAnglesMain(data_out_filename, largePertFlag, ...
                true, true, suffixStr_out) ;
            close all
        %% ---------------------------------------------------------------
        % Try to deal with problems that go with large roll perturbations
        % ---------------------------------------------------------------
        case 'extreme_roll'
            % -------------------------------------------------------------
            % load data 
            suffixStr_out = '_cleaned' ;  
            suffixStr_in  = '_cleaned' ; 
            [data_in, ~, data_out_filename, dataPath, errorFlag] = ...
                hierarchicalLoadData(pathStruct, MovNum(k), ...
                 suffixStr_out, suffixStr_in) ; 
            
             if errorFlag
                %fprintf('Could not find path for movie %d \n', MovNum(k))
                continue
             end
             % -------------------------------------------
            % back up the data folder (still not certain this works 100% of
            % the time
            fprintf('Backing up data for movie %d... \n', MovNum(k))
            [analysisPath, fn, ext] = fileparts(dataPath) ;  
            dataName = strjoin({fn, ext},'') ; 
            backupPath = fullfile(analysisPath,'backup') ; 
            if ~exist(backupPath,'dir')
                mkdir(backupPath)
            end
            copyfile(dataPath, fullfile(backupPath, dataName))
            
            % -------------------------------------------
            % perform large roll correction
            fprintf('Refining roll estimate for movie %d... \n', MovNum(k))
            [data, pseudoRollSmooth, swap_idx, check_idx] = ...
                estimatePseudoRoll(data_in, [], largePertFlag) ; 
            
            % ---------------------------
            % re-calculate angles
            fprintf('Re-calculating angles for movie %d... \n', MovNum(k))
            
            save(data_out_filename, 'data')
            data = calcAnglesMain(data_out_filename, largePertFlag, ...
                true, true, suffixStr_out) ;
            close all
        %% ---------------------------------------------------------------
        % Run script to try to fix problems with chord vector
        % ---------------------------------------------------------------
        case 'refine_chord'
            % -------------------------------------------------------------
            % load data 
            suffixStr_out = '_cleaned' ;  
            suffixStr_in  = '_cleaned' ; 
            [data_in, ~, data_out_filename, dataPath, errorFlag] = ...
                hierarchicalLoadData(pathStruct, MovNum(k), ...
                 suffixStr_out, suffixStr_in) ; 
            
             if errorFlag
                %fprintf('Could not find path for movie %d \n', MovNum(k))
                continue
             end
             % -------------------------------------------
            % back up the data folder (still not certain this works 100% of
            % the time
            fprintf('Backing up data for movie %d... \n', MovNum(k))
            [analysisPath, fn, ext] = fileparts(dataPath) ;  
            dataName = strjoin({fn, ext},'') ; 
            backupPath = fullfile(analysisPath,'backup') ; 
            if ~exist(backupPath,'dir')
                mkdir(backupPath)
            end
            copyfile(dataPath, fullfile(backupPath, dataName))
            
            % -------------------------------------------
            % perform chord refinement
            fprintf('Refining right wing chord for movie %d... \n',...
                MovNum(k))
            data = refineChordVectorSam(data_in, 'right', largePertFlag) ; 
            fprintf('Refining left wing chord for movie %d... \n',...
                MovNum(k))
            tic
            data = refineChordVectorSam(data, 'left', largePertFlag) ; 
            toc
            % ---------------------------
            % re-calculate angles
            fprintf('Re-calculating angles for movie %d... \n', MovNum(k))
            
            save(data_out_filename, 'data')
            data = calcAnglesMain(data_out_filename, largePertFlag, ...
                true, true, suffixStr_out) ;
            close all
        %% ---------------------------------------------------------------
        % Incorrect entry
        % ---------------------------------------------------------------
        otherwise
            fprintf('Invalid analysis selection for movie %d \n', MovNum(k))
            continue
    end
    %toc
    fprintf('\n Completed movie number: %d \n \n', MovNum(k))
end

% sort analysis folders into proper directories
moveEmptyExprFolders(pathToWatch)

end