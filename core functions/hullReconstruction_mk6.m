function [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
    wing1Res, wing1FrameStartInd, wing1FrameEndInd, ...
    wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
    hullReconstruction_mk6(params, CM_pos, flyBW4d, bodyBW4d, dlt, easyWandData,order)

% INPUTS
% ------
% params  - is defined in construct_params.m (script that calls this func)
% CM_pos - the center-of-mass of the body in each image. dimension of CM_pos is (3camera, Nimages, 2coordinates)
% flyBW4d
% bodyBW4d
% dlt - array of direct linear transformation coefficients. product of
%      easyWand 5
% easyWandData - includes camera principal points and focal lengths. used
%       to find transformation from world coordinates to image coordinates. also
%       a product of easyWand5
% order - vector giving the order of the cameras as given to easyWand5. 
%       order(i) = j means that camera i has index j in easyWand5. default
%       is order = [2, 1, 3]
% 
% 
% OUTPUTS
% -------
% bodyRes - Nx4 matrix for body coordinates. in the form [t x y z], where t is frame number
% bodyFrameStartInd - index of first frame in which the body is in view
% bodyFrameEndInd - index of last frame in which the body is in view
% wing1Res - Mx4 matrix for wing 1 coordinates. in the form [t x y z], where t is frame number
% wing1FrameStartInd - index of first frame in which wing 1 is in view
% wing1FrameEndInd - index of last frame in which wing 1 is in view
% wing2Res - Px4 matrix for wing 2 coordinates. in the form [t x y z], where t is frame number
% wing2FrameStartInd - index of first frame in which wing 2 is in view
% wing2FrameEndInd - index of last frame in which wing 2 is in view
% mergedWingsFlag - 1 if wings overlap, 0 otherwise
%
%
% EXAMPLE USAGE
% -------------
% [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
%     wing1Res, wing1FrameStartInd, wing1FrameEndInd, ...
%     wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
%     hullReconstruction_mk5(params, CM_pos, all_fly_bw, body_only_bw, dlt_matrix, easyWandData,[2,1,3]);

% ------------

%load(calibration_directory);
%dlt=load(dlt_directory);


if (~exist('order','var'))
    order = [ 2 1 3] ; 
    %  camera indices YZ=1 XZ=2  XY=3
    %  column indices    2    1     3
    % order(YZ) is the index of the column in dlt for YZ camera = 2 by defualt
    % order(XZ) is the index of the column in dlt for XZ camera = 1 by default
    % order(XZ) is the index of the column in dlt for XY camera = 3 by default     
end

if (isempty(order))
    order = [ 2 1 3] ; 
end

dlt_xy=dlt(:,order(3));
dlt_xz=dlt(:,order(2));
dlt_yz=dlt(:,order(1));

detectorLengthPix = params.detectorLengthPix ; % may need to change for non-square images
if numel(detectorLengthPix) == 1
    detectorLengthPix = detectorLengthPix.*[1,1] ; 
end
imageHeight = detectorLengthPix(1) ;
imageWidth  = detectorLengthPix(2) ;

%N         = params.N ;         % 512 ; % no. of voxels is N^3. % does not have to be the image size. can be larger.
volLength =  16e-3 ; % AL: was 8e-3 %params.volLength ; % N / 232 ; % 2.2 in cm

% ----------
% SET PARAMS
% ----------

YZ = params.YZ ;
XZ = params.XZ ;
XY = params.XY ;

startTrackingTime = params.startTrackingTime ; % enterTime ;
endTrackingTime   = params.endTrackingTime ; % 41 % exitTime;

Nimages           = endTrackingTime - startTrackingTime + 1 ;

% --------------------------
% MAIN LOOP - GO OVER IMAGES
% --------------------------
%indd=0;
CM_pos(:,:,2) = imageHeight - CM_pos(:,:,2) + 1 ; % convert all y coordinates to conform with easyWand requirements


[Kxz,Kyz,Kxy,Rxz,Ryz,Rxy] = find_transform2(easyWandData, dlt,order);
Nframes = endTrackingTime - startTrackingTime + 1 ;

allFlyCoords = cell(Nframes,1) ;
allBodyCoords = cell(Nframes,1) ;
allWing1Coords = cell(Nframes,1) ; 
allWing2Coords = cell(Nframes,1) ;

% prepare index offsets for 26 neighbors
ind_offsetsx = zeros(1,26,'int16') ;
ind_offsetsy = zeros(1,26,'int16') ;
ind_offsetsz = zeros(1,26,'int16') ;

c = 0;
for i3=-1:1 ;
    for j3=-1:1
        for k3=-1:1
            if ~(i3==0 && j3==0 && k3==0)
                c=c+1 ;
                ind_offsetsx(c) = i3 ;
                ind_offsetsy(c) = j3 ;
                ind_offsetsz(c) = k3 ;
            end
        end
    end
end

clear i3 j3 k3 ;

%for t=startTrackingTime:endTrackingTime
    
%    indd = indd+1 ;

firstTrackableFrame = params.firstTrackableFrame ;
NCAMS               = params.NCAMS ;
CAMERAS             = params.CAMERAS ;
voxelSize           = params.voxelSize ;

mergedWingsFlag     = false(Nframes, 1) ;
% Centers=zeros(endTrackingTime-startTrackingTime+1,3);
parfor indd= 1: (endTrackingTime-startTrackingTime+1)  %parfor
    %tic
    %disp(indd) ;
    t = indd + startTrackingTime - 1 ;
    
    %frameNum = t - startTrackingTime + 1 ;
    frameNum = t - firstTrackableFrame + 1 ; % check why this is so XXX
    
    if (indd ~= frameNum)
        error('indd ~= frameNum');
    end
    imMatFly   = false(NCAMS, imageHeight, imageWidth) ;
    imMatBody  = false(NCAMS, imageHeight, imageWidth) ;
    for c=CAMERAS
        imMatFly(c,:,:)  = getImage4D(flyBW4d, c, frameNum);
        imMatBody(c,:,:) = getImage4D(bodyBW4d, c, frameNum);
    end
    
    
    % ----------------------------------
    % MAKE SURE ALL IMAGES ARE NOT EMPTY
    % ----------------------------------
    allImagesNonEmpty = true ;
    for c=1:3
        tmp1 = squeeze(imMatFly(c,:,:)) ;
        tmp2 = squeeze(imMatBody(c,:,:)) ;
        sm1 = sum(tmp1(:)) ;
        sm2 = sum(tmp2(:)) ;
        allImagesNonEmpty = allImagesNonEmpty & (sm1>0) & (sm2>0);
    end
    
    if (~allImagesNonEmpty)
        continue ;
    end
    
    % -----------------------------------------------
    % get center of mass coordinates in each image
    % the center of mass is used to find the c.m in 3D, which is then used
    % to define the sub-volume for 3d reconstruction
    % -----------------------------------------------
    center_mass_yz = squeeze(CM_pos(YZ, indd, :))' ;
    center_mass_xz = squeeze(CM_pos(XZ, indd, :))' ;
    center_mass_xy = squeeze(CM_pos(XY, indd, :))' ; %#ok<PFBNS>
    
    % if there is no data for c.m. guess it from the image
    if (isnan(center_mass_yz))
        im = squeeze(imMatFly(YZ,:,:)) ;
        [rows, cols] = find(im) ;
        center_mass_yz(1) = mean(cols) ;
        center_mass_yz(2) = imageHeight - mean(rows) ;
    end
    
    % if there is no data for c.m. guess it from the image
    if (isnan(center_mass_xz))
        im = squeeze(imMatFly(XZ,:,:)) ;
        [rows, cols] = find(im) ;
        center_mass_xz(1) = mean(cols) ;
        center_mass_xz(2) = imageHeight - mean(rows) ;
    end
    
    % if there is no data for c.m. guess it from the image
    if (isnan(center_mass_xy))
        im = squeeze(imMatFly(XY,:,:)) ;
        [rows, cols] = find(im) ;
        center_mass_xy(1) = mean(cols) ;
        center_mass_xy(2) = imageHeight - mean(rows) ;
    end
    %clear rows cols im
    
    % calcualte for each camera separately the center of mass 2D
    % coordinates
    
    %---------------------------------------------------------------------
    % For each pair of camera/image find the camera center 3d coordinates
    % and a direction vector of the 3d line that connects
    % the camera center i with the center of mass from the image of camera i
    %-----------------------------------------------------------------------
    
    Uyz =   Ryz(1:3,1:3)' * (Kyz \ ([center_mass_yz,1]')); %#ok<PFBNS> % Ryz(1:3,1:3)' * inv(Kyz) * [center_mass_yz,1]';
    Oyz = - Ryz(1:3,1:3)' * Ryz(1:3,4);
    
    Uxz =   Rxz(1:3,1:3)' * (Kxz \ ([center_mass_xz,1]')); %#ok<PFBNS> % Rxz(1:3,1:3)' * inv(Kxz) * [center_mass_xz,1]';
    Oxz = - Rxz(1:3,1:3)' * Rxz(1:3,4);
    
    Uxy =   Rxy(1:3,1:3)' * (Kxy \ ([center_mass_xy,1]')); %#ok<PFBNS> % Rxy(1:3,1:3)' * inv(Kxy) * [center_mass_xy,1]';
    Oxy = - Rxy(1:3,1:3)' * Rxy(1:3,4);
    
    % Find G1 such that d(G1,line xz) and d(G1,line yz) are minimal
    
    A1=[-dot(Uyz,Uyz),dot(Uyz,Uxz);-dot(Uyz,Uxz),dot(Uxz,Uxz)];
    B1=[-dot(Oxz-Oyz,Uyz);-dot(Oxz-Oyz,Uxz)];
    
    lambda1=linsolve(A1,B1);
    
    G1= (Oyz+Oxz+lambda1(1)*Uyz+lambda1(2)*Uxz)/2;
    
    % Find G2 such that d(G1,line xz) and d(G1,line xy) are minimal
    
    A2=[-dot(Uxy,Uxy),dot(Uxy,Uxz);-dot(Uxy,Uxz),dot(Uxz,Uxz)];
    B2=[-dot(Oxz-Oxy,Uxy);-dot(Oxz-Oxy,Uxz)];
    
    lambda2=linsolve(A2,B2);
    
    G2= (Oxy+Oxz+lambda2(1)*Uxy+lambda2(2)*Uxz)/2;
    
    % Find G3 such that d(G3,line yz) and d(G3,line xy) are minimal
    
    A3 = [-dot(Uxy,Uxy),dot(Uxy,Uyz);-dot(Uxy,Uyz),dot(Uyz,Uyz)];
    B3 = [-dot(Oyz-Oxy,Uxy);-dot(Oyz-Oxy,Uyz)];
    
    lambda3 = linsolve(A3,B3);
    
    G3 = (Oxy+Oyz+lambda3(1)*Uxy+lambda3(2)*Uyz)/2;
    
    % Take as center the mean of the three points.
    G = [G1';G2';G3'];
    
    Ctr = mean(G); % estimated c.m. in real space.
    
%     CM_pos_mat = [center_mass_yz ; center_mass_xz; center_mass_xy] ; 
%     Ctr2 = getFlyCenterSeed(CM_pos_mat, params, dlt, voxelSize, order) ;
    %R=change_ref_frame(Rxz,Rxy);
    %Ctr=R*Ctr;
    % ---------------------------------------------
    % 3D HULL RECONSTRUCTION - FIND THE BODY VOXELS
    % ---------------------------------------------
    
    
    % shift Ctr to be the middle of some voxel.
    % given that volLength/2 / params.voxelSize is integer
    % NOTE: the points (0,0,0) in real space is the center of the camera
    % system, i.e. the center-of-mass of all the points used for
    % calibration.
    
    Ctr(1) = round(Ctr(1)/voxelSize) * voxelSize ;
    Ctr(2) = round(Ctr(2)/voxelSize) * voxelSize ;
    Ctr(3) = round(Ctr(3)/voxelSize) * voxelSize ;
    % Ctr is still in real space
    
    % find the [X0, Y0, Z0] which is the real-space coordinates of the
    % corner of the current sampling sub-volume.
    %X0=round((Ctr(1)-volLength/2)/voxelSize) * voxelSize ; % +++
    %Y0=round((Ctr(2)-volLength/2)/voxelSize) * voxelSize ;
    %Z0=round((Ctr(3)-volLength/2)/voxelSize) * voxelSize ;

    X0 = Ctr(1) - volLength/2 ; % assuming that volLength/2 is an integer number of voxelSize
    Y0 = Ctr(2) - volLength/2 ;
    Z0 = Ctr(3) - volLength/2 ;
    
    % [X0, Y0, Z0] is in the middle of some voxel

    
    % generate the real-space voxel grid vectors X, Y, Z
    % (previous code started from the middle but the X,Y,Z were the same).
    X = X0 : voxelSize : X0+volLength ;
    Y = Y0 : voxelSize : Y0+volLength ;
    Z = Z0 : voxelSize : Z0+volLength ;
    
    % find the offset of the real space volume from the real space origin
    % expressed in no. of voxels
    offset = int16([X0, Y0, Z0]/voxelSize); % +++
    
    if ( double(offset(1))*voxelSize - X0 > eps || ...
        double(offset(2))*voxelSize - Y0 > eps || ...
        double(offset(3))*voxelSize - Z0 > eps)
        error('something is wrong with voxel or offset calcuation') ;        
    end
    
    %Xind = int16( X / voxelSize - 0.5 ) ;
    %Yind = int16( Y / voxelSize - 0.5 ) ;
    %Zind = int16( Z / voxelSize - 0.5 ) ;
    
    Nx = length(X) ;
    Ny = length(Y) ;
    Nz = length(Z) ;
    
    imxyFly = squeeze(imMatFly(XY,:,:)) ; 
    imxzFly = squeeze(imMatFly(XZ,:,:)) ; 
    imyzFly = squeeze(imMatFly(YZ,:,:)) ; 
    
    imxyBody = squeeze(imMatBody(XY,:,:)) ; 
    imxzBody = squeeze(imMatBody(XZ,:,:)) ; 
    imyzBody = squeeze(imMatBody(YZ,:,:)) ; 
    
    %Q = 4 ; 
    %imyzFly(:,1:512-Q) = imyzFly(:,Q+1:512) ;
    %imyzBody(:,1:512-Q) = imyzBody(:,Q+1:512) ;
    %disp('Shifting yz up by 4 pixels') ;
    
    [imxyWing1, imxyWing2, sameMasksFlag] = segementWings(imxyFly, imxyBody) ;
    mergedWingsFlag(indd) = sameMasksFlag ;
    % define sampling volume
    %V = false(Nx,Ny,Nz);
    % tic
    % define voxel queue
    queue_q    = zeros(Nx*Ny*Nz,3, 'int16');
    queue_tail = 0 ;
    queue_head = 1 ;
    
    visited = false(Nx,Ny,Nz) ;
    goodFly = zeros(Nx*Ny*Nz,3, 'int16');
    goodFly_counter = 0 ;
    
    goodBody    = zeros(Nx*Ny*Nz,3, 'int16');
    goodBody_counter = 0 ;
    
    goodWing1    = zeros(Nx*Ny*Nz,3, 'int16');
    goodWing1_counter = 0 ;
    
    goodWing2    = zeros(Nx*Ny*Nz,3, 'int16');
    goodWing2_counter = 0 ;
    
    % find a seed voxel that is on. work with indices rather than with real
    % coordinates
    seed    = zeros(1,3) ;
    seed(1) = round(length(X)/2) ; % X( round(length(X)/2)) ;
    seed(2) = round(length(Y)/2) ; %Y( round(length(Y)/2)) ;
    seed(3) = round(length(Z)/2) ; %Z( round(length(Z)/2)) ;
    rseed = [X(seed(1)), Y(seed(2)), Z(seed(3))] ;
    if ( checkVoxel(rseed, dlt_xz, dlt_yz, dlt_xy, imxzFly, imyzFly, imxyFly, detectorLengthPix)) ;
        queue_tail = queue_tail + 1 ;
        queue_q(queue_tail,:) = seed ;        
    else % if the guess for the seed does not work
        % look in the neighbourhood
        %disp('positive seed not found in estimated c.m. position. look in the neighborhood...') ;
        neighx = round(length(X)/10) ;
        neighy = round(length(Y)/10) ;
        neighz = round(length(Z)/10) ;
        x1=max([seed(1)-neighx, 1]);
        x2=min([seed(1)+neighy, Nx]) ;
        y1=max([seed(2)-neighz, 1]);
        y2=min([seed(2)+neighx, Ny]) ;
        z1=max([seed(3)-neighy, 1]);
        z2=min([seed(3)+neighz, Nz]) ;
        found = false ;        
        for i2=x1:x2
            for j2=y1:y2 ;
                for k2=z1:z2
                    if ~found   
                        found = checkVoxel([X(i2), Y(j2), Z(k2)], dlt_xz, dlt_yz, dlt_xy, imxzFly, imyzFly, imxyFly, detectorLengthPix);
                        seed = [i2, j2, k2] ;
                        %disp('found seed!') ;
                    end
                end
            end
        end
        if ~found
            t
            figure ; imshow(imxyFly) ; title ('xy'); impixelinfo ; 
            figure ; imshow(imxzFly) ; title ('xz'); impixelinfo ; 
            figure ; imshow(imyzFly) ; title ('yz'); impixelinfo ; 
            disp('bad seed is around:') ;
            disp([X(i2), Y(j2), Z(k2)]) ;
            error('seed not found. do someting about it. maybe scan the entire volume.')
        end
        queue_tail = queue_tail + 1 ;
        queue_q(queue_tail,:) = seed ; 
    end
  
    
    while(queue_tail >= queue_head) % while queue is not empty        
        % take first item from head of queue and advance head
        curr = queue_q(queue_head,:) ;
        queue_head = queue_head + 1 ;

        visited(curr(1), curr(2), curr(3)) = true ;
        
        rcurr = [X(curr(1)), Y(curr(2)), Z(curr(3))] ;
        
        
        [isFly, isBody, isWing1, isWing2] = ...
            checkVoxelEverything(rcurr, dlt_xz, dlt_yz, dlt_xy, ...
                                imxzFly, imyzFly, imxyFly, ...
                                imxzBody, imyzBody, imxyBody, ...
                                imxyWing1, imxyWing2, detectorLengthPix); 
        
        if (sameMasksFlag && isWing1 && ~isWing2)
            isWing2 = true ;
        end
        
        if (sameMasksFlag && isWing2 && ~isWing1)
            isWing1 = true ;
        end
        
%          if (indd==102)
%              if (rcurr(1)<-0.010 && rcurr(2)>=0.01800 && rcurr(3)>-0.0053)
%                  if (isWing1 || isWing2)
%                      WWW = 'ro' ;
%                  else
%                      WWW = 'k.' ; 
%                  end
%                 figure(99) ; hold on ; plot3(rcurr(1), rcurr(2), rcurr(3),WWW) ;
%                 pause(0.001) ;
%              end
%              
%          end
        
        %if (checkVoxel(rcurr, dlt_xz, dlt_yz, dlt_xy, imxz, imyz, imxy, detectorLengthPix))
        if (isFly)
            % add current voxel to the good list
            goodFly_counter = goodFly_counter + 1 ;
            goodFly(goodFly_counter,:) = curr ;
            
            if (isBody)
                goodBody_counter = goodBody_counter + 1 ;
                goodBody(goodBody_counter,:) = curr ;
            end
            
            if (isWing1)
                goodWing1_counter = goodWing1_counter + 1 ;
                goodWing1(goodWing1_counter,:) = curr ;
                %disp(goodWing1_counter)
                %{
                if curr(1) == 0 && curr(2) == 0 && curr(3) == 0
                    keyboard ;
                end
                if goodWing1_counter == 4672
                    figure(99) ;
                    hold on
                    plot3(goodWing1(1:(goodWing1_counter-1), 1),goodWing1(1:(goodWing1_counter-1), 2),goodWing1(1:(goodWing1_counter-1), 3), 'ko','MarkerSize',8)
                    disp(goodWing1_counter)
                    keyboard ;
                end
                %}
            end
            
            if (isWing2)
                goodWing2_counter = goodWing2_counter + 1 ;
                goodWing2(goodWing2_counter,:) = curr ;
            end
            
            
            % check neighbors and add the unvisited ones to the q
            indx = curr(1) + ind_offsetsx ;
            indy = curr(2) + ind_offsetsy ;
            indz = curr(3) + ind_offsetsz ;
            
            % trim the voxels outside the box
            ii = indx>0 & indx<=Nx & indy>0 & indy<=Ny & indz>0 & indz<=Nz ;
            indx = indx(ii) ;
            indy = indy(ii) ;
            indz = indz(ii) ;
            for p=1:length(indx)
                if ~visited(indx(p), indy(p), indz(p))
                    % if not visited, add to q
                    queue_tail = queue_tail + 1 ;
                    queue_q(queue_tail,:) = [indx(p), indy(p), indz(p)] ;
                    visited(indx(p), indy(p), indz(p)) = true ;
                end
            end
        end
        
        % check if it's on
        % if so, add it to the good list and add its neighbors to the q
 
    end
  
    currFlyCoords   = goodFly(1:goodFly_counter,:) ;
    currBodyCoords  = goodBody(1:goodBody_counter,:) ;
    currWing1Coords = goodWing1(1:goodWing1_counter,:) ;
    currWing2Coords = goodWing2(1:goodWing2_counter,:) ;
    
    if (0)
        figure(78) ; %#ok<UNRCH>
        subplot(1,3,1) ; imshow(imyzFly) ; hold on ; plot(center_mass_yz(1), imageHeight - center_mass_yz(2),'ro') ; title('yz') ;
        subplot(1,3,2) ; imshow(imxyFly) ; hold on ; plot(center_mass_xy(1), imageHeight - center_mass_xy(2),'ro') ; title('xy') ;
        subplot(1,3,3) ; imshow(imxzFly) ; hold on ; plot(center_mass_xz(1), imageHeight - center_mass_xz(2),'ro') ; title('xz') ;
    end
     

    %{
    count = 0;
    
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                %{
                pixelxz    = dlt_inverse(dlt_xz, [X(i),Y(j),Z(k)] );
                pixelxz(1) = round(pixelxz(1));
                pixelxz(2) = round(512-pixelxz(2)+1);
                
                if ~(pixelxz(1)>0 && pixelxz(1)<detectorLengthPix && pixelxz(2)>0 && pixelxz(2)<detectorLengthPix ...
                        && imxz(pixelxz(2),pixelxz(1)))
                    continue ;
                end
                
                pixelyz    = dlt_inverse(dlt_yz, [X(i),Y(j),Z(k)]);
                pixelyz(1) = round(pixelyz(1));
                pixelyz(2) = round(512-pixelyz(2)+1);
                
                if ~(pixelyz(1)>0 && pixelyz(1)<detectorLengthPix && pixelyz(2)>0 && pixelyz(2)<detectorLengthPix ...
                        && imyz(pixelyz(2),pixelyz(1)))
                    continue ;
                end
                
                pixelxy    = dlt_inverse(dlt_xy, [X(i), Y(j), Z(k)]);
                pixelxy(1) = round(pixelxy(1));
                pixelxy(2) = round(512-pixelxy(2)+1);
                
                if (pixelxy(1)>0 && pixelxy(1)<detectorLengthPix && pixelxy(2)>0 && pixelxy(2)<detectorLengthPix ...
                        && imxy(pixelxy(2),pixelxy(1)))
                    V(i,j,k) = true ;
                    count = count + 1 ;
                end
                %}
                if (checkVoxel([X(i), Y(j), Z(k)]))
                    V(i,j,k) = true ;
                    count = count + 1 ;
                end
            end
        end
    end
    
    currCoords     = zeros(count,3, 'int16');
    currCoordsReal = zeros(count,3, 'double');
    ind=0;
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                if (V(i,j,k))
                    ind=ind+1;
                    currCoords(ind,:)     = [Xind(i), Yind(j), Zind(k)];
                    currCoordsReal(ind,:) = [X(i),       Y(j),    Z(k)];
                end
            end
        end
    end
    
    %}
    
    %allFlyCoords{indd} = currFlyCoords ;
    
    % here the body and wing voxel coordinates are stored. they are offset
    % by the "offset" vector, such that the voxels are measured with
    % respect to the lab frame of reference, for which (0,0,0) is the center of the
    % overall filming volume.
    allBodyCoords{indd}  = currBodyCoords  + repmat(offset,size(currBodyCoords, 1), 1); % +++
    allWing1Coords{indd} = currWing1Coords + repmat(offset,size(currWing1Coords,1), 1);
    allWing2Coords{indd} = currWing2Coords + repmat(offset,size(currWing2Coords,1), 1);
    %if isempty(allWing1Coords{indd})
    %    allWing1Coords{indd} = [nan nan nan] ;
    %end
    %if isempty(allWing2Coords{indd})
    %    allWing2Coords{indd} = [nan nan nan] ;
    %end
    
    currFlyCoordsReal   = zeros(goodFly_counter,  3, 'double');
    currBodyCoordsReal  = zeros(goodBody_counter, 3, 'double') ;
    currWing1CoordsReal = zeros(goodWing1_counter,3, 'double') ;
    currWing2CoordsReal = zeros(goodWing2_counter,3, 'double') ;
    
    for p=1:goodFly_counter
        currFlyCoordsReal(p,:) = [X(currFlyCoords(p,1)), Y(currFlyCoords(p,2)), Z(currFlyCoords(p,3))] ;
    end
    
    for p=1:goodBody_counter
        currBodyCoordsReal(p,:) = [X(currBodyCoords(p,1)), Y(currBodyCoords(p,2)), Z(currBodyCoords(p,3))] ;
    end
    
    for p=1:goodWing1_counter
        currWing1CoordsReal(p,:) = [X(currWing1Coords(p,1)), Y(currWing1Coords(p,2)), Z(currWing1Coords(p,3))] ;
    end
    
    for p=1:goodWing2_counter
        currWing2CoordsReal(p,:) = [X(currWing2Coords(p,1)), Y(currWing2Coords(p,2)), Z(currWing2Coords(p,3))] ;
    end
    
    %toc
    %keyboard
    if (0)
        figure(77) ;
        set(gcf,'position',[212   400   921   565]) ; 
        clf ; hold on ;
        %plot3(currFlyCoordsReal(:,1), currFlyCoordsReal(:,2), currFlyCoordsReal(:,3),'bo','markersize',2) ; grid on ; box on ;
        plot3(currBodyCoordsReal(:,1), currBodyCoordsReal(:,2), currBodyCoordsReal(:,3),'go','markersize',2) ; 
        plot3(currWing1CoordsReal(:,1), currWing1CoordsReal(:,2), currWing1CoordsReal(:,3),'ro','markersize',2) ; 
        plot3(currWing2CoordsReal(:,1), currWing2CoordsReal(:,2), currWing2CoordsReal(:,3),'bo','markersize',2) ; 
        grid on ; box on ;
        plot3(Ctr(1), Ctr(2), Ctr(3),'ro','markerfacecolor','r') ;
        plot3(X(1), Y(1), Z(1),'k^') ;
        plot3(X(end), Y(end), Z(end),'k^') ;
        hold on ; axis equal ; axis tight ;
        title(['t = ' num2str(t) ]) ;
        view(-25,34) ;
        keyboard ; 
    end
    
    if(0)
        ww = [-1 1 -1 1] * 40  ; %#ok<UNRCH>
        figure(77); clf ;
        subplot(2,2,1) ;
        imshow(imyzFly) ; title('yz') ; 
        xc = center_mass_yz(1) ;
        yc = imageHeight - center_mass_yz(2) ;        
        hold on ; plot(xc, yc,'ro') ; title('yz') ; axis([xc xc yc yc]+ww) ;
      
        subplot(2,2,2) ;
        imshow(imxzFly) ; title('xz') ; 
        xc = center_mass_xz(1) ;
        yc = imageHeight - center_mass_xz(2) ;        
        hold on ; plot(xc, yc,'ro') ; title('xz') ; axis([xc xc yc yc]+ww) ;
        
        subplot(2,2,3) ;        
        imshow(imxyFly) ; title('xy') ; 
        xc = center_mass_xy(1) ;
        yc = imageHeight - center_mass_xy(2) ;        
        hold on ; plot(xc, yc,'ro') ; title('xy') ; axis([xc xc yc yc]+ww) ;
        
        
        subplot(2,2,4) ;
        plot3(currFlyCoordsReal(:,1), currFlyCoordsReal(:,2), currFlyCoordsReal(:,3),'bo','markersize',2) ; grid on ; box on ;
        hold on ;
        plot3(Ctr(1), Ctr(2), Ctr(3),'ro','markerfacecolor','r') ;
        plot3(X(1), Y(1), Z(1),'k^') ;
        plot3(X(end), Y(end), Z(end),'k^') ;
        hold on ; axis equal
        title(['t = ' num2str(t) ]) ;
        view(-71,8) ;
        %keyboard;
    end 
    %toc
end

% calculate size to allocate for bodyRes
Sbody = zeros(Nimages,1,'double') ;
Swing1 = zeros(Nimages,1,'double') ;
Swing2 = zeros(Nimages,1,'double') ;

for n=1:Nimages
    Sbody(n) = size(allBodyCoords{n},1) ;
    
    if size(allWing1Coords{n},1) >= 1        %changed by SW, 7/2/15
        Swing1(n) = size(allWing1Coords{n},1) ;
    else
        Swing1(n) = 1 ;
    end
    
    if size(allWing2Coords{n},1) >= 1       %changed by SW, 7/2/15
        Swing2(n) = size(allWing2Coords{n},1) ;
    else
        Swing2(n) = 1 ;
    end
end

bodyRes = zeros(sum(Sbody), 4,'int16');
wing1Res = zeros(sum(Swing1), 4,'int16');
wing2Res = zeros(sum(Swing2), 4,'int16');

nextIndBody = 1 ;
nextIndWing1 = 1 ;
nextIndWing2 = 1 ;

for t=startTrackingTime:endTrackingTime
    n = t - startTrackingTime + 1 ;
    
    bodyRes(nextIndBody:nextIndBody+Sbody(n)-1,:) = [ int16(ones(Sbody(n),1)*t) allBodyCoords{n} ];
    nextIndBody = nextIndBody + Sbody(n) ;
    
    if Swing1(n) > 1 %changed by SW, 7/2/15
        wing1Res(nextIndWing1:nextIndWing1+Swing1(n)-1,:) = [ int16(ones(Swing1(n),1)*t) allWing1Coords{n} ];
    else
        wing1Res(nextIndWing1:nextIndWing1+Swing1(n)-1,:) = [ int16(ones(Swing1(n),1)*t) mean([mean(allWing1Coords{n-1}); mean(allWing1Coords{n+1})]) ] ;
    end
    nextIndWing1 = nextIndWing1 + Swing1(n) ;
    
    if Swing2(n) > 1 %changed by SW, 7/2/15
        wing2Res(nextIndWing2:nextIndWing2+Swing2(n)-1,:) = [ int16(ones(Swing2(n),1)*t) allWing2Coords{n} ];
    else
        wing2Res(nextIndWing2:nextIndWing2+Swing2(n)-1,:) = [ int16(ones(Swing2(n),1)*t) mean([mean(allWing2Coords{n-1}); mean(allWing2Coords{n+1})]) ] ;
    end
    nextIndWing2 = nextIndWing2 + Swing2(n) ;
end

%{
[ bodyRes, bodyFrameStartInd, bodyFrameEndInd ...
    wing1Res, wing1FrameStartInd, wing1FrameEndInd ...
    wing2Res, wing2FrameStartInd, wing2FrameEndInd ] = .

%}

df = diff(bodyRes(:,1)) ;
bodyFrameStartInd = [1 ; find(df==1)+1] ;
bodyFrameEndInd   = [bodyFrameStartInd(2:end)-1 ; size(bodyRes,1)] ;

df = diff(wing1Res(:,1)) ;
wing1FrameStartInd = [1 ; find(df==1)+1] ;
wing1FrameEndInd   = [wing1FrameStartInd(2:end)-1 ; size(wing1Res,1)] ;

df = diff(wing2Res(:,1)) ;
wing2FrameStartInd = [1 ; find(df==1)+1] ;
wing2FrameEndInd   = [wing2FrameStartInd(2:end)-1 ; size(wing2Res,1)] ;

%change the reference frame 

%framesRange   = [startTrackingTime endTrackingTime] ;

end % function hullReconstruction_mk5(...)

% ------------------------------------------------------------------------

function bool = checkVoxel(r, dlt_xz, dlt_yz, dlt_xy, imxz, imyz, imxy, detectorLengthPix) % was a "nested function", not it's not (due to parfor)

bool = false ; % default value
%r1=R'*r';
%r=r1';
pixelxz    = dlt_inverse(dlt_xz, r );
pixelxz(1) = round(pixelxz(1));
pixelxz(2) = round(detectorLengthPix(1)-pixelxz(2)+1);

if ~(pixelxz(1)>0 && pixelxz(1)<detectorLengthPix(2) && pixelxz(2)>0 && ...
        pixelxz(2)<detectorLengthPix(1) && imxz(pixelxz(2),pixelxz(1)))
    return
end


pixelyz    = dlt_inverse(dlt_yz, r);
pixelyz(1) = round(pixelyz(1));
pixelyz(2) = round(detectorLengthPix(1)-pixelyz(2)+1);

if ~(pixelyz(1)>0 && pixelyz(1)<detectorLengthPix(2) && pixelyz(2)>0 && ...
        pixelyz(2)<detectorLengthPix(1) && imyz(pixelyz(2),pixelyz(1)))
    return
end

pixelxy    = dlt_inverse(dlt_xy, r);
pixelxy(1) = round(pixelxy(1));
pixelxy(2) = round(detectorLengthPix(1)-pixelxy(2)+1);

if (pixelxy(1)>0 && pixelxy(1)<detectorLengthPix(2) && pixelxy(2)>0 && ...
        pixelxy(2)<detectorLengthPix(1) && imxy(pixelxy(2),pixelxy(1)))
    bool = true ;
end
end % function checkVoxel(r)



function [isFly, isBody, isWing1, isWing2] = checkVoxelEverything(r, dlt_xz, dlt_yz, dlt_xy, ...
    imxzFly, imyzFly, imxyFly, ...
    imxzBody, imyzBody, imxyBody, ...
    imxyWing1, imxyWing2, detectorLengthPix) 


isFly   = false ; % default values
isBody  = false ;
isWing1 = false ;
isWing2 = false ;

% first check if the voxel belongs to the "entire fly", same as in
% checkVoxel. If yes, check if it belongs to the body or wings (or none).
%r1=R'*r';
%r=r1';
pixelxz    = dlt_inverse(dlt_xz, r );
pixelxz(1) = round(pixelxz(1));
pixelxz(2) = round(detectorLengthPix(1)-pixelxz(2)+1);

if ~(pixelxz(1)>0 && pixelxz(1)<detectorLengthPix(2) && ...
        pixelxz(2)>0 && pixelxz(2)<detectorLengthPix(1) ...
        && imxzFly(pixelxz(2),pixelxz(1)))
    return
end

pixelyz    = dlt_inverse(dlt_yz, r);
pixelyz(1) = round(pixelyz(1));
pixelyz(2) = round(detectorLengthPix(1)-pixelyz(2)+1);

if ~(pixelyz(1)>0 && pixelyz(1)<detectorLengthPix(2) && ...
        pixelyz(2)>0 && pixelyz(2)<detectorLengthPix(1) ...
        && imyzFly(pixelyz(2),pixelyz(1)))
    return
end

pixelxy    = dlt_inverse(dlt_xy, r);
pixelxy(1) = round(pixelxy(1));
pixelxy(2) = round(detectorLengthPix(1)-pixelxy(2)+1);

if (pixelxy(1)>0 && pixelxy(1)<detectorLengthPix(2) && ...
        pixelxy(2)>0 && pixelxy(2)<detectorLengthPix(1) ...
        && imxyFly(pixelxy(2),pixelxy(1)))
    isFly = true ;
else
    return
end

% if we got to this line, then the voxel belongs to the entire fly. now
% check if it belongs to the body, wing1, wing2, or none of them.

% check for body - iff voxel belongs to imxyBody, imxzBody and imyzBody

if ( imxyBody(pixelxy(2),pixelxy(1)) && ...
        imyzBody(pixelyz(2),pixelyz(1)) && ...
        imxzBody(pixelxz(2),pixelxz(1))) 
    isBody = true ;
    return
end

% if we got to this line, then the voxel is not body.

flag = imyzBody(pixelyz(2),pixelyz(1)) & imxzBody(pixelxz(2),pixelxz(1)) ;

if (imxyWing1(pixelxy(2),pixelxy(1)) && ~flag)
    isWing1 = true ;
    return
end


if (imxyWing2(pixelxy(2),pixelxy(1)) && ~flag)
    isWing2 = true ;
    return
end



end % function checkVoxelEverything(r)

% -----------------------------------------------------------------------
function [Kxz,Kyz,Kxy,Rxz,Ryz,Rxy] = find_transform2(easyWandData, dlt,order)

%load(calibration_directory)
%dlt=load(coefs_directory);
%order is a 1x3 vector  such that order(1)= index of camera YZ for calibration
%                                 order(2)= index of camera XZ for calibration
%                                 order(3)= index of camera XY for calibration  

dltxy=dlt(:,order(3));
dltxz=dlt(:,order(2));
dltyz=dlt(:,order(1));

pp=easyWandData.ppts;
f=easyWandData.focalLengths';

f=f(:,order);
pp=pp(:,[2*order(1)-1,2*order(1),2*order(2)-1,2*order(2),2*order(3)-1,2*order(3)]);

Kyz=[f(1),0,pp(1);0,f(1),pp(2);0,0,1];
Kxz=[f(2),0,pp(3);0,f(2),pp(4);0,0,1];
Kxy=[f(3),0,pp(5);0,f(3),pp(6);0,0,1];

Axz=[dltxz(1),dltxz(2),dltxz(3),dltxz(4);dltxz(5),dltxz(6),dltxz(7),dltxz(8);dltxz(9),dltxz(10),dltxz(11),1];
Ayz=[dltyz(1),dltyz(2),dltyz(3),dltyz(4);dltyz(5),dltyz(6),dltyz(7),dltyz(8);dltyz(9),dltyz(10),dltyz(11),1];
Axy=[dltxy(1),dltxy(2),dltxy(3),dltxy(4);dltxy(5),dltxy(6),dltxy(7),dltxy(8);dltxy(9),dltxy(10),dltxy(11),1];

Rxz= Kxz \ Axz ; % inv(Kxz)*Axz;
Ryz= Kyz \ Ayz ; % inv(Kyz)*Ayz;
Rxy= Kxy \ Axy ; % inv(Kxy)*Axy;

Rxz=Rxz/max([norm(Rxz(1:3,1)),norm(Rxz(1:3,2)),norm(Rxz(1:3,3))]);
Ryz=Ryz/max([norm(Ryz(1:3,1)),norm(Ryz(1:3,2)),norm(Ryz(1:3,3))]);
Rxy=Rxy/max([norm(Rxy(1:3,1)),norm(Rxy(1:3,2)),norm(Rxy(1:3,3))]);


end


% ------------------------------------------------------------------------
function [imWing1, imWing2, sameMasksFlag] = segementWings(imFly, imBody,...
    WING_AREA_THRESH)

if ~exist('WING_AREA_THRESH', 'var')
    WING_AREA_THRESH = 12 ;
end

sameMasksFlag = false ;
oneWingEmptyFlag = false ;

thickbody = bwmorph(imBody,'thicken',1) ;
dif       = imFly - thickbody ;
dif       = bwmorph(dif,'clean') ; % remove isolated pixels
% find connected componenets and sort by size

CC  = bwconncomp(dif);
Ncc = length(CC.PixelIdxList) ;
    
svec = zeros(Ncc,1) ; % vector containing the size of each connected components
for j=1:Ncc
    svec(j) = length(CC.PixelIdxList{j}) ;
end

[svec, idx] = sort(svec,'descend') ;
numObjects  = length(idx) ;
wingsIdx    = zeros(2,1) ;

if (numObjects>=2)
    for j=1:2
        if (svec(j)>=WING_AREA_THRESH)
            wingsIdx(j) = idx(j) ;
        else
            wingsIdx(j) = -1 ;
        end
    end
elseif (numObjects==1)
    wingsIdx = [ 1 ; 1 ] ;
    sameMasksFlag = true ;
else
    disp('No objects... hmm...') ;
    keyboard ;
end

imWing1 = false(size(imFly)) ;
imWing2 = false(size(imFly)) ;

% keep only large objects
if (wingsIdx(1)>0)
    imWing1(CC.PixelIdxList{wingsIdx(1)}) = true ;
    %imWing1(objData(wingsIdx(1)).PixelIdxList) = true ; % old code
end

if (wingsIdx(2)>0)
    imWing2(CC.PixelIdxList{wingsIdx(2)}) = true ;
    %imWing2(objData(wingsIdx(2)).PixelIdxList) = true ; % old code
end

% if at least one wing is empty, use dif instead
% later on the two wings will be resolved using clustering in
% hullAnalysis_mk?
if (wingsIdx(1)==-1 || wingsIdx(2)==-1)
    
    imWing1 = false(size(dif)) ;
    
    [label, numObjects] = bwlabel(dif, 4);
    if (numObjects>1)
        objData = regionprops(label, 'area','PixelIdxList') ;
        [~, idx] = max([objData.Area]) ;
        imWing1(objData(idx).PixelIdxList) = true ;
    else
        imWing1 = dif ;
    end
    imWing2 = imWing1 ;
    sameMasksFlag = true ;
end

end

% ------------------------------------------------------------------------
%% generate camera ray passing through pixel in 3D space using DLT
% generate points O and V which are the coordinates for the camera center
% and the center of an input pixel (respectively) in 3D space
%
% NB: pix_val (1x2 array of pixel coordinates) should be in easyWand form,
% i.e. do the detectorLengthPix - pix_val(2) + 1
%
function [O, U] = myCameraRay(pix_val, dlt, camNum, order)
if ~exist('order','var')
    order = [2, 1, 3] ;
end

dlt_mat = reshape([dlt(:, order(camNum)) ; 1],[4,3])' ; % get dlt in matrix
U  = pinv(dlt_mat)*[pix_val,1]';                        % pixel coordinate
O  = null(dlt_mat) ;                                    % camera center
U = U(1:3)./U(4) ; O = O(1:3)./O(4) ;                   % normalize
end

% ------------------------------------------------------------------------
%% shortest distance between two skew lines (wikipedia)
% v1, v2 are arbitray points on line v
% u1, u2 are arbitray points on line u
function dist = skewLineDist(v1, v2, u1, u2)

% write lines as vectors
p1 = v1 ;
d1 = (v2-v1)./norm(v2-v1) ;

p2 = u1 ;
d2 = (u2-u1)./norm(u2-u1) ;

% find new line that is perpendicular to both lines
n = cross(d1, d2) ;

% calculate distance
dist = abs(dot(n, (p2 - p1))) ;

end

% -------------------------------------------------------------------------
%% shortest line segment joining two skew lines (wikipedia)
% v1, v2 are arbitray points on line v
% u1, u2 are arbitray points on line u
function [c1, c2] = shortestConnLine(v1, v2, u1, u2)

% write lines as vectors
p1 = v1 ;
d1 = (v2-v1)./norm(v2-v1) ;

p2 = u1 ;
d2 = (u2-u1)./norm(u2-u1) ;

% find new line that is perpendicular to both lines
n = cross(d1, d2) ;

% plane formed by translations of u along n contains p2 and is perp to n2:
n2 = cross(d2, n) ;
n1 = cross(d1, cross(d2, d1)) ;

% find points
c1 = p1 + (dot((p2-p1), n2)/dot(d1,n2)).*d1 ;
c2 = p2 + (dot((p1-p2), n1)/dot(d2,n1)).*d2 ;

end

% -------------------------------------------------------------------------
%% get fly center
function Ctr = getFlyCenterSeed(CM_pos_mat, params, dlt, voxelSize, order)
if ~exist('order','var') || isempty(order)
    order = [2, 1, 3] ;
end
%---------------------------------------------------------------------
% For each pair of camera/image find the camera center 3d coordinates
% and a direction vector of the 3d line that connects
% the camera center i with the center of mass from the image of camera i
%-----------------------------------------------------------------------
YZ = params.YZ ; XZ = params.XZ ; XY = params.XY ; % camera indices

% data storage
U_mat = nan(3) ; % Ndim (xyz) x Ncams. vectors for pixel center of mass
O_mat = nan(3) ; % Ndim (xyz) x Ncams. vectors for camera center
for c = 1:size(CM_pos_mat,1)
    [O_mat(:,c), U_mat(:,c)] = myCameraRay(CM_pos_mat(c,:), dlt, c, order) ;
end

% loop through the various ray comparions we want to make
G = nan(3) ; % Ndim (xyz) x Ncams. vectors for points that minimize distance between rays
comparison_list = [ XZ, YZ ; XZ, XY ; YZ, XY] ;
for i = 1:size(comparison_list,1)
    comp_idx = comparison_list(i,:) ;
    % Find G1 such that d(G1,line xz) and d(G1,line yz) are minimal
    [line_pt1, line_pt2] = shortestConnLine(U_mat(:,comp_idx(1)), ...
        O_mat(:,comp_idx(1)), U_mat(:,comp_idx(2)), O_mat(:,comp_idx(2))) ;
    G(i,:) = (line_pt1' + line_pt2')./2 ;
end

% Take as center the mean of the three points.
Ctr = mean(G); % estimated c.m. in real space.

% shift Ctr to be the middle of some voxel.
% given that volLength/2 / params.voxelSize is integer
% NOTE: the points (0,0,0) in real space is the center of the camera
% system, i.e. the center-of-mass of all the points used for
% calibration.

Ctr(1) = round(Ctr(1)/voxelSize) * voxelSize ;
Ctr(2) = round(Ctr(2)/voxelSize) * voxelSize ;
Ctr(3) = round(Ctr(3)/voxelSize) * voxelSize ;

% check with plots?
if (0)
    figure ;
    hold on
    color_cell = {'r', 'g', 'b', 'm', 'c', 'y'} ;
    for j = 1:size(O_mat, 2)
        plot3([O_mat(1,j), U_mat(1,j)], [O_mat(2,j), U_mat(2,j)], ...
            [O_mat(3,j), U_mat(3,j)],'-','LineWidth',2, ...
            'Color', color_cell{j})
    end
    plot3(Ctr(1), Ctr(2), Ctr(3), 'ko','MarkerFaceColor','k')
    axis equal
    grid on
    box on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    set(gca,'xlim', Ctr(1) + 50*voxelSize*[-1, 1])
    set(gca,'ylim', Ctr(2) + 50*voxelSize*[-1, 1])
    set(gca,'zlim', Ctr(3) + 50*voxelSize*[-1, 1])
end

end