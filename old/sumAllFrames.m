
metaData = getCinMetaData(cinFilename) ;

if(~exist('tin','var'))
    tin = metaData.firstImage ;
    tout = metaData.lastImage ;
else
    if (isempty(tin))
        tin = metaData.firstImage ;
        tout = metaData.lastImage ;
    end
end

cindata = myOpenCinFile(cinFilename) ;

firstIm =  myReadCinImage(cindata, tin) ;
summedImage = zeros(size(firstIm)) ;
se = strel('disk',6) ; 

c = 0 ;
for t=tin:tout
    
    c=c+1 ;
    
    %im1 = ReadCineFileImage(cinFilename, t, false);
    im1 = myReadCinImage(cindata, t) ;
    
    im2 = imsubtract(bg, im1) ;
    im3 = imopen(im2,se) ; 
    summedImage = summedImage + double(im2) ; 
end

summedImageNormalized = summedImage./max(summedImage(:)) ;
figure ; imshow(summedImageNormalized,[]) ;

Iadj = imadjust(summedImageNormalized) ;
figure ; imshow(Iadj,[])
Iadj2 = imopen(Iadj,se) ;
figure ; imshow(Iadj2) 

level = graythresh(Iadj2) ;
bw = im2bw(Iadj2,level) ;

bw2 = bwareaopen(bw,40) ;
figure ; imshow(bw)