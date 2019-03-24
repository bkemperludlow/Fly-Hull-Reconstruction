%--------------------------------------------------------------------------
% script to change indices on cine/xml files from flight experiments, since
% the numbering is sometimes offset due to multiple files being saved in
% quick succession and/or the data acqusition code renumbering the files in
% an incorrect way
%
% this differs from "reNumberMovieFiles.m" in that it tries to globally 
% renumber all movies as opposed to just renaming some. this script also
% takes longer to run
%
% %formerly called "change_nb_v3.m"
%--------------------------------------------------------------------------
vidPath = 'D:\Box Sync\VNC Sensory Lines\17_09122018\to_be_sorted\' ;
pathCurr = pwd ; 
cd (vidPath)

%cinCell = cell(1,3) ; %xy = 1, xz = 2, yz = 3
%xmlCell = cell(1,3) ;
datenum_tol = (1/60)*(1/24) ; %should correspond to 1 minute, in units of days 

xyCinDir = dir('xy*.cin') ; xzCinDir = dir('xz*.cin') ; yzCinDir = dir('yz*.cin') ;
xyXMLDir = dir('xy*.xml') ; xzXMLDir = dir('xz*.xml') ; yzXMLDir = dir('yz*.xml') ;

numXYCin = length(xyCinDir) ;
numXZCin = length(xzCinDir) ;
numYZCin = length(yzCinDir) ;
numCinArray = [numXYCin, numXZCin, numYZCin] ; %xy = 1, xz = 2, yz = 3
numCinAll = sum(numCinArray) ; 

xyCinDatenum = [xyCinDir(:).datenum]' ; xyXMLDatenum = [xyXMLDir(:).datenum]' ; 
xzCinDatenum = [xzCinDir(:).datenum]' ; xzXMLDatenum = [xzXMLDir(:).datenum]' ; 
yzCinDatenum = [yzCinDir(:).datenum]' ; yzXMLDatenum = [yzXMLDir(:).datenum]' ; 

datenumCinAll = [ xyCinDatenum , ones(size(xyCinDatenum)) ; ...
    xzCinDatenum , 2*ones(size(xzCinDatenum)) ; ...
    yzCinDatenum , 3*ones(size(yzCinDatenum))] ; 

% datenumXMLAll = [ xyXMLDatenum , ones(size(xyXMLDatenum)) ; ...
%    xzXMLDatenum , 2*ones(size(xzXMLDatenum)) ; ...
%    yzXMLDatenum , 3*ones(size(yzXMLDatenum))] ; 

[~, cinSortInd] = sort(datenumCinAll(:,1)) ;
datenumCinAll_sorted = datenumCinAll(cinSortInd,:) ; 

% [~, xmlSortInd] = sort(datenumXMLAll(:,1)) ;
% datenumXMLAll_sorted = datenumXMLAll(xmlSortInd,:) ; 

cinDirAll = [xyCinDir ; xzCinDir ; yzCinDir ] ; 
cinDirAll_sorted = cinDirAll(cinSortInd) ; 
xmlDirAll = [xyXMLDir ; xzXMLDir ; yzXMLDir] ; 
xmlDirAll_sorted = xmlDirAll(xmlSortInd) ; 

alreadySorted = [] ; 

for i = 1:numCinAll 
    
    if sum(alreadySorted == i) > 0 
        continue ; 
    end
    
    datenum_curr = datenumCinAll_sorted(i,1) ; 
    within_tol_ind = (abs(datenumCinAll_sorted(:,1) - datenum_curr) < datenum_tol) ;  
    
    if length(unique(datenumCinAll_sorted(within_tol_ind,2))) == 3
        xy_ind = find((datenumCinAll_sorted(:,2) == 1) & within_tol_ind , 1, 'first') ;
        xz_ind = find((datenumCinAll_sorted(:,2) == 2) & within_tol_ind, 1, 'first') ;
        yz_ind = find((datenumCinAll_sorted(:,2) == 3) & within_tol_ind, 1, 'first') ;
        
        alreadySorted = [alreadySorted ; xy_ind ; xz_ind ; yz_ind] ; 
    else
        continue ; 
    end
    
    
end

mkdir([vidPath 'sorted\'])
cd([vidPath 'sorted\']) 

cc = 0 ;

for j = 1:3:length(alreadySorted)
    for m = 0:2
        [~, f, extstr] = fileparts(cinDirAll_sorted(alreadySorted(j+m)).name);
        copyfile(strcat(vidPath, f, extstr), strcat(f(1:3), num2str(cc,'%03u'), '.cin'));
        
%         [~, g, extstr2] = fileparts(xmlDirAll_sorted(alreadySorted(j+m)).name);
%         copyfile(strcat(rootName, g, extstr2), strcat(g(1:3), num2str(cc,'%03u'), '.xml'));
    end
    cc = cc + 1 ; 
end

%--------------------------------------------------------------------------
% return to original folder
cd(pathCurr)

%{
if (numXYCin ~= numXZCin) || (numXYCin ~= numYZCin) || (numYZCin ~= numXZCin)
    disp('Different number of movies for the different cameras')
    [numCin, camRef] = min(numCinArray) ;
    dateNumCell = cell(1, 3) ;
    goodIndArray = zeros(numCin,3) ; %xy = 1, xz = 2, yz = 3
    
    dateNumCell{1} = xyCinDir(:).datenum ; 
    dateNumCell{2} = xzCinDir(:).datenum ; 
    dateNumCell{3} = yzCinDir(:).datenum ; 
    
    for j = 1:3
        if j == camRef
            goodIndArray(:,j) = 1:numCin ;  
            continue ;
        else
            dateDist = pdist2(dateNumCell{camRef},dateNumCell{j}) ;
            minDateDist = min(dateDist) ;
            [~, sortInd] = sort(minDateDist,'ascend') ;
            
            goodIndArray(:,j) = sortInd(1:numCin) ;
        end
    end
    diffCinNumFlag = true ;
else
    numCin = numXYCin ;
    %goodIndArray = repmat([1:numCin]',1,3) ; %xy = 1, xz = 2, yz = 3
    diffCinNumFlag = false ;
    [~,xyCinSortInd] = sort([xyCinDir(:).datenum]) ; xyCinDir = xyCinDir(xyCinSortInd) ; 
    [~,xySortInd] = sort([xyXMLDir(:).datenum]) ; xyXMLDir = xyXMLDir(xyCinSortInd) ; 
    [~,xzCinSortInd] = sort([xzCinDir(:).datenum]) ; xzCinDir = xzCinDir(xyCinSortInd) ; 
    [~,xzXMLSortInd] = sort([xzXMLDir(:).datenum]) ; xzXMLDir = xzXMLDir(xyCinSortInd) ; 
    [~,yzCinSortInd] = sort([yzCinDir(:).datenum]) ; yzCinDir = yzCinDir(xyCinSortInd) ; 
    [~,yzXMLSortInd] = sort([yzXMLDir(:).datenum]) ; yzXMLDir = yzXMLDir(xyCinSortInd) ; 
    
end    
    
cinDir = [xyCinDir(goodIndArray(:,1)), xzCinDir(goodIndArray(:,2)), yzCinDir(goodIndArray(:,3)) ] ;
xmlDir = [xyXMLDir(goodIndArray(:,1)), xzXMLDir(goodIndArray(:,2)), yzXMLDir(goodIndArray(:,3)) ] ;

mkdir([rootName 'sorted\'])
cd([rootName 'sorted\']) 

cc = 0 ;

for id = 1:numCin
    
    if cc<10
        zstr = '00' ;
    elseif cc<100
        zstr = '0' ;
    else
        zstr = '' ;
    end
    
    for m = 1:3
        [~, f, extstr] = fileparts(cinDir(id,m).name);
        copyfile(strcat(rootName, f, extstr), strcat(f(1:3), zstr, num2str(goodIndArray(cc,m),3), '.cin'));
    
        [~, g, extstr2] = fileparts(xmlDir(id,m).name);
        copyfile(strcat(rootName, g, extstr2), strcat(g(1:3), zstr, num2str(goodIndArray(cc,m),3), '.xml'));
    end
    
    cc = cc + 1;
    %str=['ren ' filePrefix '_' zstr1 int2str(i) fileExtension ' ' filePrefix '_' zstr2 int2str(i+deltaInd) fileExtension];
    
    %dos (str);
end


%str=['ren xy_0' int2str(i) '.cin xy_0' int2str(i-1) '.cin'];
%str=['ren ' filePrefix '_00' int2str(i) fileExtension ' ' filePrefix '_00' int2str(i+deltaInd) fileExtension];
%}