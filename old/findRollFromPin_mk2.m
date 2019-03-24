function [rollAngle] = findRollFromPin_mk2(tip,data)

%Estimates the roll angle in each frame of a movie using the
%tip vector

%Developed by Luca and Sam
%
%Input
%-----
%data - the output of hullAnalysis_mk3. Notably contains the coordinates of
%the body voxels of the fly, as well as the body axis

%tip- the output of find_pin_mk3. Contains the tip coordinates of the pin
%for each frame of the movie.

%Output
%rollAngle- the roll angle in each frame



%% First initialize array and load relevant aspects of 'data'
Nframes = data.Nimages ;
Ntimes = size(tip(:,4),1);
bodyCM=data.bodyCM;
AHat=data.AHat;
rollAngle = zeros(Ntimes,1);


%FIND THE TIP POSITION IN AN UNROLLED FRAME

%FrameRef=input('Enter the reference frame number frame where roll angle is 0\n');
FrameRef= 351;
refAHat = AHat(FrameRef,:);
% refCM = bodyCM(1,:);

idx=find(tip(:,4)==FrameRef);

refTipVect = tip(idx,1:3)-bodyCM(FrameRef,:);
rollAngle(idx) = 0;
pitchAngleRef=(180/pi)*asin(AHat(FrameRef,3));
yawAngleRef=(180/pi)*atan2(AHat(FrameRef,2),AHat(FrameRef,1));

% pitchAngleRef=pitchAngles(FrameRef);
% yawAngleRef=yawAngles(FrameRef);

%Unyaw then unpitch the tip to have his coordinates
% Same thing for body axis

refTipVect=RotatePoint(refTipVect,[0 0 0], [0 0 1], -yawAngleRef);
refTipVect=RotatePoint(refTipVect,[0 0 0],[0 1 0],pitchAngleRef); %check that it is the correct axis.
refTipVect=refTipVect/norm(refTipVect);

refAHat=RotatePoint(refAHat,[0 0 0], [0 0 1],-yawAngleRef);
refAHat=RotatePoint(refAHat,[0 0 0],[0 1 0],pitchAngleRef);
%check that refAhat is aligned with x.



%% Find the roll angle in each frame
for i=1:Nframes
    
    if isempty(find(tip(:,4)==i))
        continue;
    end
    
    if i==FrameRef
        continue;
    end
    
    %if isnan(tip(i,1))
    %    continue;
    %end
    
    idx=find(tip(:,4)==i);
    
    currTipVect=tip(idx,1:3)-bodyCM(i,:);
    %     currTipVect=tip(i,:);
    yawAngle=(180/pi)*atan2(AHat(i,2),AHat(i,1));
    pitchAngle=(180/pi)*asin(AHat(i,3));
    %     yawAngle=yawAngles(i);
    %     pitchAngle=pitchAngles(i);
    currTipVect=RotatePoint(currTipVect,[0 0 0], [0 0 1],-yawAngle);
    %     compTipVect=RotatePoint(refTipVect,[0 0 0],[0 0 1],-pitchAngle);
    currTipVect=RotatePoint(currTipVect,[0 0 0],[0 1 0],pitchAngle);
    
    %normalization of both vectors
    currTipVect=currTipVect/norm(currTipVect);
    %     compTipVect=compTipVect/norm(compTipVect);
    projCurrTipVect = currTipVect-dot(currTipVect,[1 0 0])*[1 0 0];
    projrefTipVect = refTipVect-dot(refTipVect,[1 0 0])*[1 0 0];
    
    projCurrTipVect=projCurrTipVect/norm(projCurrTipVect);
    projrefTipVect=projrefTipVect/norm(projrefTipVect);
    
    %estimate rollAngle
    
    %
    %     projCurrTipVect = currTipVect-dot(currTipVect,[0 0 1])*[0 0 1];
    %     projCompTipVect = compTipVect-dot(compTipVect,[0 0 1])*[0 0 1];
    
         checkvect=cross(projrefTipVect,projCurrTipVect);
    
         if checkvect(1)>0
             rollAngle(idx)= acosd(dot(projCurrTipVect,projrefTipVect));
    
         else
             rollAngle(idx)= -acosd(dot(projCurrTipVect,projrefTipVect));
         end
    
end

rollAngle = [rollAngle, tip(:,4)];
figure;plot(rollAngle(:,2),rollAngle(:,1));
end