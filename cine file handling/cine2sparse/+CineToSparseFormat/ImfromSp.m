function [Im] = ImfromSp(fr,metaData,frames)

Im=zeros(metaData.frameSize);
IndIm=sub2ind(metaData.frameSize,frames(fr).indIm(:,1),frames(fr).indIm(:,2));
Im(IndIm)=frames(fr).indIm(:,3);

end

