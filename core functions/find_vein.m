function [rightVein,leftVein ] = find_vein( data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rightSpan=data.rightSpanHats;
leftSpan=data.leftSpanHats;
rightChord=data.rightChordHats;
leftChord=data.leftChordHats;
Nimages=size(rightSpan,1);
rightVein=zeros(Nimages,3);
leftVein=zeros(Nimages,3);

for i=1:Nimages

rightAxis=cross(rightSpan(i,:),rightChord(i,:));
leftAxis=cross(leftSpan(i,:),leftChord(i,:));

rightVein(i,:)=RotatePoint(rightSpan(i,:),[0 0 0],rightAxis,12);
leftVein(i,:)=RotatePoint(leftSpan(i,:),[0 0 0],leftAxis,12);
end

end

