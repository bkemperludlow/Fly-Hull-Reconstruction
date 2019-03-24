rootName = 'F:\Sam\Janelia Flies\08_MB258C UAS-TNT-E\' ; %include backslash!
cd (rootName)

Nfolders = 0 ;
existFlag = 1 ;
folderNames = {rootName} ;

%checks for folders within a larger experiment folder of the form "run 2"
%etc.
while existFlag > 0 
    folderStr = ['run ' num2str(Nfolders+1,1)] ;
    if exist(folderStr) ~= 7
        existFlag = existFlag - 1 ;
    end
    Nfolders = Nfolders + 1 ;
end

%loads directories for xy, yz, and xz for both cine and xml files into cell
%arrays from original folder
cinCell = cell(Nfolders,3) ; %xy = 1, xz = 2, yz = 3
xmlCell = cell(Nfolders,3) ;

cinCell{1,1} = dir('xy*.cin') ; cinCell{1,2} = dir('xz*.cin') ; cinCell{1,3} = dir('yz*.cin') ;
xmlCell{1,1} = dir('xy*.xml') ; xmlCell{1,2} = dir('xz*.xml') ; xmlCell{1,3} = dir('yz*.xml') ;
Nmovies = length(cinCell{1,1}) ;

%repeats above process for any other folders (e.g. run 2)
if Nfolders > 2
    for i = 2:Nfolders
        newFolderName = [rootName 'run ' num2str(i-1) '\'] ;
        folderNames = [folderNames, {newFolderName}] ; 
        cd(newFolderName)
        cinCell{i,1} = dir('xy*.cin') ; cinCell{i,2} =  dir('xz*.cin') ; cinCell{i,3} = dir('yz*.cin') ;
        xmlCell{i,1} = dir('xy*.xml') ; xmlCell{i,2} =  dir('xz*.xml') ; xmlCell{i,3} = dir('yz*.xml') ;
        Nmovies = Nmovies + length(cinCell{i,1}) ;
    end
end    

% make new directory to store sorted movies--this way old ones don't get
% overwritten
cd (rootName)
mkdir([rootName 'sorted\'])
cd([rootName 'sorted\']) 

fileDate = zeros(Nmovies,1) ;
c = 1 ;
for j = 1:Nfolders
    Ntemp = length(cinCell{j,1}) ;
    for id = 1:Ntemp 
        fileDate(c) = cinCell{j,1}(id).datenum ;
        c = c + 1 ; 
    end
end

[~,I] = sort(fileDate) ;
c2 = 1 ;

for k = 1:Nfolders
    pathstr = folderNames(k) ;
    Ntemp = length(cinCell{k,1}) ;
    for id = 1:Ntemp
        for m = 1:3
            if I(c2) < 10
                zstr = '00' ;
            elseif I(c2) < 100
                zstr = '0' ;
            else
                zstr = '' ;
            end
            
            [~, f, extstr] = fileparts(cinCell{k,m}(id).name);
            copyfile(strcat(pathstr{1}, f, extstr), strcat(f(1:3), zstr, num2str(I(c2),3), '.cin'));
            
            [~, g, extstr2] = fileparts(xmlCell{k,m}(id).name);
            copyfile(strcat(pathstr{1}, g, extstr2), strcat(g(1:3), zstr, num2str(I(c2),3), '.xml'));
        end
        c2 = c2 + 1 ;
    end
    
end
   





