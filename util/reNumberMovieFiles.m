%--------------------------------------------------------------------------
% script to change indices on cine/xml files from flight experiments, since
% the numbering is sometimes offset due to multiple files being saved in
% quick succession and/or the data acqusition code renumbering the files in
% an incorrect way
%
% this differs from "reNumberMovieFilesAll.m" in that it only alters movie
% names by adding/subtracting a fixed number from specified files, rather
% than try to globally renumber all movies
%
% Formerly called "change_nb.m"
%--------------------------------------------------------------------------
vidPath = [] ; % 'D:\Box Sync Old\Opto Silencing\46_23102020\' ; 
pathCurr = pwd ;
cd(vidPath)

% which camera files to rename
filePrefixArray =  {'xy', 'xz', 'yz'} ; % {'xy', 'xz', 'yz'} ; % {'xy', 'xz', 'yz'} ;
% which file types to rename
fileExtensionArray = {'.cine', '.xml'} ;

% keyboard

% define the start, end, and shift number for the movie files you'd like to
% rename
startInd = 0 ; 
endInd = 0 ; 
deltaInd = +1 ;
% used to prevent overwriting
backwardsFlag = true ; 

for n = 1:length(filePrefixArray)
    %set camera
    filePrefix = filePrefixArray{n} ;
    for m = 1:length(fileExtensionArray)
        %set file type
        fileExtension = fileExtensionArray{m} ;
        for j= startInd:endInd
            % switch indices if we're going in reverse numerical order
            if backwardsFlag
                i = endInd + startInd - j ; 
            else
                i = j ;
            end
            
            % generate string for dos command to run that renames files
            str=['ren ' filePrefix '_' num2str(i,'%03d') fileExtension ' ' ...
                filePrefix '_' num2str(i+deltaInd,'%03d') fileExtension];
            
            % run dos command
            dos (str);
        end
    end
end

% return to original folder
cd(pathCurr) 
%str=['ren xy_0' int2str(i) '.cin xy_0' int2str(i-1) '.cin'];
%str=['ren ' filePrefix '_00' int2str(i) fileExtension ' ' filePrefix '_00' int2str(i+deltaInd) fileExtension];