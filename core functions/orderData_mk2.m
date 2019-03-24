
%{

Original version of HRMT saved the following variables

% body
stroke(:,2) = xs(:);
stroke(:,3) = ys(:);
stroke(:,4) = zs(:);
stroke(:,5) = (180/pi)*psis(:);
stroke(:,6) = (180/pi)*betas(:);
stroke(:,7) = (180/pi)*rhos(:);

% right wing
stroke(:,8) = x1s(:);
stroke(:,9) = y1s(:);
stroke(:,10) = z1s(:);
stroke(:,11) = (180/pi)*phi1s(:);
stroke(:,12) = (180/pi)*theta1s(:);
stroke(:,13) = (180/pi)*eta1s(:);

% left wing
stroke(:,14) = x2s(:);
stroke(:,15) = y2s(:);
stroke(:,16) = z2s(:);
stroke(:,17) = (180/pi)*phi2s(:);
stroke(:,18) = (180/pi)*theta2s(:);
stroke(:,19) = (180/pi)*eta2s(:);

% user corrections
flags = zeros(numImages,19);
adjust = zeros(numImages,19);
swap = zeros(numImages,1);
flip1 = zeros(numImages,1);
flip2 = zeros(numImages,1);

MY CODE
-------

hullReconstruction results were stored in a res-type matrix, 
(e.g. load res1to532.mat)


initial results from hullAnalysis.m were saved in a file according to:

save (outputFilename, ...
        'RESIDX', 'stroke', 'rollHats', 'span1Hats', 'chord1Hats', 'chord1Hat2s', ...
        'phi1Hats', 'span2Hats', 'chord2Hats', 'chord2Hat2s', 'phi2Hats', 'bodyInd', 'rightWingInd', ...
        'psiHats', 'normRolls', 'leftWingInd', 'Nimages','params', 'errorLog') ;

(e.g. load hullAnalysis_1_to_532_.mat)

the corrections made by the user in the script originalGUIscript2.m are:

1. the original varialbe used by HRMT is the "adjust" matrix as well as 
   "flag", "swap" and "flip" matrices (see above).
   My GUI script RE-calculates the body and wing angles according to the 
   user's VECTOR adjustments, but this is done in the lab frames. To clarify,
   the adjustments in "adjust" take into account non-commutative rotations,
   
2. My GUI script also stores the VECTORS which are the result of the user's
   corrections. These used to re-calc the angles (see 1 above) and ALSO
   to calc the angles in the BODY frame of reference (done in the temporary
   function plotHullAnalysisResults.m). By the way it is much 
   simpler to do the lab frame calculations using vectors, since we can
   rotate them easily. 
   The five vectors are:

    chord1Hats_user
    chord2Hats_user 
    span1Hats_user
    span2Hats_user
    rollHats_user

(saved in file e.g. load userCorrections_1_to_532.mat)

What is missing?
----------------
1. Storing all data in a convenient structure
   a. stored data should probably include user corrections in the "final"
      data rather than keeping original data and user corrections separately.

   b. data to save (store angles or only vectors from which to calc angles?)
      
      res
      RESIDX

      Body center of mass: xs ys zs
      Right wing center of mass: x1s x2s x3s
      Left wing center of mass: x1s x2s x3s

      5 vectors that include user corrections: 
      chord1Hats, chord2Hats, span1Hats, span2Hats, rollHats

      pin parameters (returned by findPinPositionManual.m):
      rollAngles, lineParams, validPoints, timeIndices, userInput

      misc: params, errorLog (still need to define this array)

2. calculate the observables with respect to either
   a. Lab frame of reference
   b. Body frame of reference.

3. plot the observables

4. figure out how to calibrate pin orientation with respect to the body



%}

%{
usage for example:
mov34data = orderData_mk2(res, RESIDX, stroke, 0, ...
chord1Hats, chord2Hats, span1Hats, span2Hats, rollHats, ...
chord1AltHats, chord2AltHats, ...
-1, -1, -1, -1, params, errorLog) ;
%}

% improvedAHat: [532x3 double]
% improvedAHatMidPoints: [532x3 double]

function orderedData = orderData_mk2 (res, RESIDX, stroke, adjust, ...
    chord1Hats, chord2Hats, span1Hats, span2Hats, rollHats, ...
    chord1AltHats, chord2AltHats, ...
    bodyRollAngles, pinLineParams, pinValidPoints, pinUserInput, ...
    params, errorLog)

Nimages = size(stroke,1) ;

orderedData.Nimages = Nimages ;

stroke = stroke + adjust ;

% voxel identities
orderedData.res = res ;
orderedData.RESIDX = RESIDX ;

% constants
orderedData.bodyInd      = 1 ;
orderedData.rightWingInd = 2 ;
orderedData.leftWingInd  = 3 ;

% center of mass
orderedData.bodyCM      = stroke(:,2:4  ) ;
orderedData.rightWingCM = stroke(:,8:10 ) ;
orderedData.leftWingCM  = stroke(:,14:16) ;

% five vectors 
orderedData.rightChordHats = chord1Hats ;
orderedData.leftChordHats  = chord2Hats ;
orderedData.rightSpanHats  = span1Hats ;
orderedData.leftSpanHats   = span2Hats ;
orderedData.AHat           = rollHats ;

orderedData.chord1AltHats = chord1AltHats ;
orderedData.chord2AltHats = chord2AltHats ;

orderedData.rollVectors    = zeros(Nimages, 3) ;

% pin params
orderedData.pin.bodyRollAngles = bodyRollAngles ;
orderedData.pin.pinLineParams  = pinLineParams ; 
orderedData.pin.pinValidPoints = pinValidPoints ;
orderedData.pin.pinUserInput   = pinUserInput ;

% misc
orderedData.params   = params ;
orderedData.errorLog = errorLog ;

orderedData.rhoTimes = [] ;


end

