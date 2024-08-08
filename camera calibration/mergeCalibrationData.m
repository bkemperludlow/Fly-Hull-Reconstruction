% file type
ext = '.cine' ; % '.cin'

% cam info
camNames = {'xy', 'xz', 'yz'} ; 
camOrder = [2, 3, 1] ; % NB: this is just to do with order within this script -- should probably change to make more consistent
Ncams = length(camNames) ; 

% path to calibration movies
calibPath = pwd ; 
if ~contains(calibPath,'calibration')
   fprintf('Warning: calibration path seems incorrect \n')
   keyboard
end

% params
rad1 = 17 ; 
rad2 = 50 ; % 45 ; 

% find calibration videos in directory
calibDir = dir(fullfile(calibPath, ['*',ext])) ; 
camCineIdx = cell(Ncams,1) ; 
for c = 1:Ncams
    camCineIdx{c} = find(arrayfun(@(y) contains(y.name, camNames{c},...
        'IgnoreCase',1), calibDir)) ; 
end

% get wand points for each round of videos
Ncines = length(camCineIdx{1}) ; 
Mcell = cell(Ncines,1) ; % initialize wand point storage

% loop over videos
for n = 1:Ncines
   Mcell{n} = getWandPoints(rad1, rad2, calibDir(camCineIdx{1}(n)).name, ...
       calibDir(camCineIdx{2}(n)).name, calibDir(camCineIdx{3}(n)).name,...
       sprintf('wandPoints%d.csv',n)) ; 
   
   fprintf('Completed calibration movie %d/%d \n', n, Ncines)
end


% close parallel pool
delete(gcp) 

% compile wandpoints into one array
M = vertcat(Mcell{:}) ;  
dlmwrite('wandPoints.csv', M) ; % write to file

% ----------------------------------------------------
% generate camera profile file
focalLengthEst = 7000 ;  % assume a fixed focal length estimage

% get image size and principle points (loop over cameras)
imSizes = zeros(3,2) ; 
for c = 1:Ncams
    ind = camOrder(c) ; 
    cine_fn = calibDir(camCineIdx{ind}(1)).name ;
    md = getCinMetaData(cine_fn) ;
    
    imSizes(c, :) = [md.height, md.width] ; 
end
principalPoints = round(imSizes./2) ;

% generate file
fid = fopen('fake_camera_profile.txt','wt');
for c = 1:Ncams
    % ind = camOrder(c) ;
    entry_curr = fprintf(fid, '%d %d %d %d %d %d 1 0 0 0 0 0\n', c, ...
        focalLengthEst, imSizes(c,2), imSizes(c,1), ...
        principalPoints(c,2), principalPoints(c,1)) ;
end
fclose(fid);

% run easyWand
easyWand5

