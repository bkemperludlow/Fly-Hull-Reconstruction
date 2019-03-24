rootName = 'G:\Janelia Flies\08_MB258C UAS-TNT-E\' ;
cd (rootName)

%cinCell = cell(1,3) ; %xy = 1, xz = 2, yz = 3
%xmlCell = cell(1,3) ;

xyCinDir = dir('xy*.cin') ; xzCinDir = dir('xz*.cin') ; yzCinDir = dir('yz*.cin') ;
xyXMLDir = dir('xy*.xml') ; xzXMLDir = dir('xz*.xml') ; yzXMLDir = dir('yz*.xml') ;

numXYCin = length(xyCinDir) ;
numXZCin = length(xzCinDir) ;
numYZCin = length(yzCinDir) ;
numCinArray = [numXYCin, numXZCin, numYZCin] ; %xy = 1, xz = 2, yz = 3

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
    diffCinNumFLag = true ;
else
    numCin = numXYCin ;
    goodIndArray = repmat([1:numCin]',1,3) ; %xy = 1, xz = 2, yz = 3
    diffCinNumFLag = false ;
end    
    
cinDir = [xyCinDir(goodIndArray(:,1)), xzCinDir(goodIndArray(:,2)), yzCinDir(goodIndArray(:,3)) ] ;
xmlDir = [xyXMLDir(goodIndArray(:,1)), xzXMLDir(goodIndArray(:,2)), yzXMLDir(goodIndArray(:,3)) ] ;

mkdir([rootName 'sorted\'])
cd([rootName 'sorted\']) 

cc = 1 ;

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